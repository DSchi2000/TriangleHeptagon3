// 9-variable rational phase-space generator over GF(p), p = 2^31-1
//
// Reads Roots.txt — one expression per line; all must be QR mod p.
// Expressions are arbitrary rational functions of the variables (not just
// polynomials): +, -, *, /, ^, parentheses, and division by a variable or
// any sub-expression are all supported, e.g.
//     (s12*s23 - s34^2) / (c1 + c2)^2
//
// Unlike gen_6pt, there is no Gram constraint to solve here: the
// parameterization {s12, s23, s34, s123, s234, s56, c1, c2, c3} already
// satisfies G6 = 0 identically, so every sampled point already lies on the
// constraint surface — no aux variable, no quadratic solve.
//
// Two sampling modes:
//   full sample (default): every variable drawn independently and uniformly
//   --slice <a> <b>:        1-parameter line s = a + t*b, t random in GF(p)
//
// Adding/removing roots: just edit Roots.txt — no recompilation needed.

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <omp.h>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

// ── Prime ─────────────────────────────────────────────────────────────────────
static constexpr int64_t P = 2147483647LL; // 2^31 - 1

// ── Variable table ────────────────────────────────────────────────────────────
static const char* const VAR_NAMES[] = {
    "s12","s23","s34","s123","s234","s56","c1","c2","c3"
};
static constexpr int N_VARS = 9;

// ── Modular arithmetic ────────────────────────────────────────────────────────
static inline int64_t addmod(int64_t a, int64_t b) {
    int64_t r = a + b; return r >= P ? r - P : r;
}
static inline int64_t submod(int64_t a, int64_t b) {
    return a >= b ? a - b : a + P - b;
}
static inline int64_t mulmod(int64_t a, int64_t b) {
    return (int64_t)((__int128)a * b % P);
}
static int64_t powmod(int64_t b, int64_t e) {
    int64_t r = 1; b %= P; if (b < 0) b += P;
    for (; e > 0; e >>= 1, b = mulmod(b, b)) if (e & 1) r = mulmod(r, b);
    return r;
}
static inline int64_t invmod(int64_t a)  { return powmod(a, P - 2); }
static inline int64_t sqrtmod(int64_t a) { return a == 0 ? 0 : powmod(a, (P + 1) / 4); }

// Jacobi symbol (a / P)
static int jacobi(int64_t a, int64_t n) {
    a %= n; if (a < 0) a += n;
    int j = 1;
    while (a) {
        while (!(a & 1)) { a >>= 1; int64_t r = n & 7; if (r==3||r==5) j=-j; }
        std::swap(a, n);
        if ((a & 3) == 3 && (n & 3) == 3) j = -j;
        a %= n;
    }
    return n == 1 ? j : 0;
}
static inline bool is_square(int64_t a) { return jacobi(a, P) >= 0; }

// ── Expression AST (rational functions of the variables) ─────────────────────
struct Node {
    enum Kind { NUM, VAR, ADD, SUB, MUL, DIV, NEG, POW } kind;
    int64_t value = 0;   // NUM
    int     var_idx = -1;// VAR
    int     exp = 0;     // POW (may be negative)
    std::unique_ptr<Node> a, b;
};
using NodePtr = std::unique_ptr<Node>;

static NodePtr make_num(int64_t v) {
    auto n = std::make_unique<Node>(); n->kind = Node::NUM;
    n->value = ((v % P) + P) % P; return n;
}
static NodePtr make_var(int idx) {
    auto n = std::make_unique<Node>(); n->kind = Node::VAR; n->var_idx = idx; return n;
}
static NodePtr make_bin(Node::Kind k, NodePtr a, NodePtr b) {
    auto n = std::make_unique<Node>(); n->kind = k; n->a = std::move(a); n->b = std::move(b); return n;
}
static NodePtr make_neg(NodePtr a) {
    auto n = std::make_unique<Node>(); n->kind = Node::NEG; n->a = std::move(a); return n;
}
static NodePtr make_pow(NodePtr a, int e) {
    auto n = std::make_unique<Node>(); n->kind = Node::POW; n->a = std::move(a); n->exp = e; return n;
}

// Returns false iff a division by zero (mod P) was hit — the expression is
// undefined at this point, so the caller should treat the sample as invalid.
static bool eval_node(const Node* n, const int64_t* vars, int64_t& out) {
    switch (n->kind) {
        case Node::NUM: out = n->value; return true;
        case Node::VAR: out = vars[n->var_idx]; return true;
        case Node::NEG: {
            int64_t x;
            if (!eval_node(n->a.get(), vars, x)) return false;
            out = submod(0, x); return true;
        }
        case Node::ADD: {
            int64_t x, y;
            if (!eval_node(n->a.get(), vars, x) || !eval_node(n->b.get(), vars, y)) return false;
            out = addmod(x, y); return true;
        }
        case Node::SUB: {
            int64_t x, y;
            if (!eval_node(n->a.get(), vars, x) || !eval_node(n->b.get(), vars, y)) return false;
            out = submod(x, y); return true;
        }
        case Node::MUL: {
            int64_t x, y;
            if (!eval_node(n->a.get(), vars, x) || !eval_node(n->b.get(), vars, y)) return false;
            out = mulmod(x, y); return true;
        }
        case Node::DIV: {
            int64_t x, y;
            if (!eval_node(n->a.get(), vars, x) || !eval_node(n->b.get(), vars, y)) return false;
            if (y == 0) return false;
            out = mulmod(x, invmod(y)); return true;
        }
        case Node::POW: {
            int64_t x;
            if (!eval_node(n->a.get(), vars, x)) return false;
            if (n->exp >= 0) { out = powmod(x, n->exp); return true; }
            int64_t p = powmod(x, -n->exp);
            if (p == 0) return false;
            out = invmod(p); return true;
        }
    }
    return false; // unreachable
}

static int node_count(const Node* n) {
    if (!n) return 0;
    return 1 + node_count(n->a.get()) + node_count(n->b.get());
}

// ── Parser: expr := term (('+'|'-') term)*
//             term := factor (('*'|'/') factor)*
//             factor := ('+'|'-')* primary ('^' ['-'] INT)?
//             primary := INT | IDENT | '(' expr ')'
static int var_index(const char* s, int len) {
    for (int i = 0; i < N_VARS; i++) {
        const char* n = VAR_NAMES[i];
        if ((int)strlen(n) == len && memcmp(n, s, len) == 0) return i;
    }
    throw std::runtime_error("Unknown variable: " + std::string(s, len));
}

struct Lex {
    const char* p;
    const char* end;
    void skip()  { while (p < end && isspace((unsigned char)*p)) ++p; }
    char peek()  { skip(); return p < end ? *p  : 0; }
    bool eat(char c) { if (peek()==c) { ++p; return true; } return false; }
    bool at_end()   { skip(); return p >= end; }

    int64_t read_int() {
        skip(); int64_t v = 0;
        if (p >= end || !isdigit((unsigned char)*p))
            throw std::runtime_error("Expected a number");
        while (p < end && isdigit((unsigned char)*p)) v = v * 10 + (*p++ - '0');
        return v;
    }
    std::pair<const char*, int> read_ident() {
        skip(); const char* s = p;
        while (p < end && (isalnum((unsigned char)*p) || *p == '_')) ++p;
        return {s, (int)(p - s)};
    }
};

static NodePtr parse_expr(Lex& lx);

static NodePtr parse_primary(Lex& lx) {
    if (lx.eat('(')) {
        auto e = parse_expr(lx);
        if (!lx.eat(')')) throw std::runtime_error("Expected ')'");
        return e;
    }
    char c = lx.peek();
    if (isdigit((unsigned char)c)) return make_num(lx.read_int());
    if (isalpha((unsigned char)c)) {
        auto [s, len] = lx.read_ident();
        return make_var(var_index(s, len));
    }
    throw std::runtime_error("Unexpected character in expression");
}

static NodePtr parse_factor(Lex& lx) {
    bool neg = false;
    while (true) {
        if (lx.eat('-')) { neg = !neg; continue; }
        if (lx.eat('+')) { continue; }
        break;
    }
    auto base = parse_primary(lx);
    if (lx.eat('^')) {
        bool eneg = lx.eat('-');
        int e = (int)lx.read_int();
        base = make_pow(std::move(base), eneg ? -e : e);
    }
    if (neg) base = make_neg(std::move(base));
    return base;
}

static NodePtr parse_term(Lex& lx) {
    auto left = parse_factor(lx);
    while (true) {
        if (lx.eat('*')) left = make_bin(Node::MUL, std::move(left), parse_factor(lx));
        else if (lx.eat('/')) left = make_bin(Node::DIV, std::move(left), parse_factor(lx));
        else break;
    }
    return left;
}

static NodePtr parse_expr(Lex& lx) {
    auto left = parse_term(lx);
    while (true) {
        if (lx.eat('+')) left = make_bin(Node::ADD, std::move(left), parse_term(lx));
        else if (lx.eat('-')) left = make_bin(Node::SUB, std::move(left), parse_term(lx));
        else break;
    }
    return left;
}

static NodePtr parse_line(const std::string& line) {
    const char* p = line.c_str();
    size_t n = line.size();
    // Strip an optional leading "<digits>\t" equation-number prefix (common
    // when pasting expressions out of a CAS).
    const char* q = p;
    while (q - p < 12 && q < p + n && isdigit((unsigned char)*q)) ++q;
    if (q > p && q < p + n && *q == '\t') { p = q + 1; n = line.size() - (p - line.c_str()); }

    Lex lx{p, p + n};
    auto e = parse_expr(lx);
    if (!lx.at_end())
        throw std::runtime_error("Unexpected trailing text in: " + line);
    return e;
}

static std::vector<NodePtr> load_roots(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("Cannot open: " + path);
    std::vector<NodePtr> roots;
    std::string line;
    int lineno = 0;
    while (std::getline(f, line)) {
        ++lineno;
        if (line.empty() || line[0] == '#') continue;
        try {
            roots.push_back(parse_line(line));
        } catch (const std::exception& e) {
            throw std::runtime_error(path + ":" + std::to_string(lineno) + ": " + e.what());
        }
    }
    return roots;
}

// ── Load vector from file ─────────────────────────────────────────────────────
static std::vector<int64_t> load_vec(const std::string& path, int n) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("Cannot open: " + path);
    std::vector<int64_t> out(n);
    for (int i = 0; i < n; i++) {
        if (!(f >> out[i]))
            throw std::runtime_error(path + ": expected " + std::to_string(n) + " values");
        out[i] = ((out[i] % P) + P) % P;
    }
    return out;
}

// ── Usage ─────────────────────────────────────────────────────────────────────
static void usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " [options]\n"
        << "  --roots   <path>    Roots.txt              (default: Roots.txt)\n"
        << "  --npoints <N>       points to generate      (default: 100)\n"
        << "  --seed    <S>       RNG seed                 (default: 42)\n"
        << "  --threads <N>       number of CPU cores      (default: 1)\n"
        << "  --out     <file>    also save to file\n"
        << "  --sqrt              also output square roots of each expression\n"
        << "  --slice <a> <b>     slice mode: s = a + t*b for random t,\n"
        << "                      vectors must cover all " << N_VARS << " variables:\n"
        << "                      s12 s23 s34 s123 s234 s56 c1 c2 c3\n"
        << "  N                   positional: equivalent to --npoints N\n"
        << "\n"
        << "No Gram constraint is solved here — the parameterization already\n"
        << "satisfies it, so with no --slice every variable is just sampled\n"
        << "independently and uniformly.\n";
}

// ── Main ─────────────────────────────────────────────────────────────────────
int main(int argc, char* argv[]) {
    std::string roots_path   = "Roots.txt";
    int64_t     target       = 100;
    uint64_t    seed         = std::random_device{}();
    int         n_threads    = 1;
    bool        do_sqrts     = false;
    std::string out_path;
    std::string slice_a_path, slice_b_path;

    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        if      (a == "--roots"   && i+1 < argc) roots_path   = argv[++i];
        else if (a == "--out"     && i+1 < argc) out_path     = argv[++i];
        else if (a == "--npoints" && i+1 < argc) target       = std::stoll(argv[++i]);
        else if (a == "--seed"    && i+1 < argc) seed         = std::stoull(argv[++i]);
        else if (a == "--threads" && i+1 < argc) n_threads    = std::stoi(argv[++i]);
        else if (a == "--slice"   && i+2 < argc) { slice_a_path = argv[++i]; slice_b_path = argv[++i]; }
        else if (a == "--sqrt")  do_sqrts = true;
        else if (a == "--help")  { usage(argv[0]); return 0; }
        else if (isdigit(a[0])) target = std::stoll(a);
    }

    bool slice_mode = !slice_a_path.empty();
    std::vector<int64_t> vec_a, vec_b;
    if (slice_mode) {
        try {
            vec_a = load_vec(slice_a_path, N_VARS);
            vec_b = load_vec(slice_b_path, N_VARS);
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << "\n"; return 1;
        }
        std::cerr << "Slice mode: s = a + t*b  (t random in GF(p))\n";
    } else {
        std::cerr << "Full-sample mode: all " << N_VARS << " variables drawn independently\n";
    }

    std::vector<NodePtr> roots;
    try {
        roots = load_roots(roots_path);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n"; return 1;
    }
    std::cerr << "Seed: " << seed << "\n";
    std::cerr << "Loaded " << roots.size() << " root conditions\n";
    std::cerr << "Threads: " << n_threads << "\n";

    // Cheapest-first check order (by AST node count) — doesn't own the nodes.
    std::vector<const Node*> check_order;
    for (const auto& r : roots) check_order.push_back(r.get());
    std::sort(check_order.begin(), check_order.end(),
              [](const Node* x, const Node* y) { return node_count(x) < node_count(y); });

    // Optional output file
    std::ofstream outfile;
    if (!out_path.empty()) {
        outfile.open(out_path);
        if (!outfile) { std::cerr << "Error: cannot open " << out_path << "\n"; return 1; }
        std::cerr << "Writing to: " << out_path << "\n";
    }

    std::mutex output_mutex;
    auto writeln = [&](const std::string& line) {
        std::lock_guard<std::mutex> lock(output_mutex);
        std::cout << line << "\n";
        if (outfile) outfile << line << "\n";
    };

    // Header
    {
        std::string hdr = "#";
        for (int i = 0; i < N_VARS; i++) { hdr += " "; hdr += VAR_NAMES[i]; }
        if (slice_mode) hdr += " t";
        if (do_sqrts)
            for (size_t i = 0; i < roots.size(); i++)
                hdr += " sqrt" + std::to_string(i + 1);
        writeln(hdr);
    }

    std::atomic<int64_t> found{0}, trials{0};

    #pragma omp parallel num_threads(n_threads)
    {
        std::mt19937_64 local_rng(seed + (uint64_t)omp_get_thread_num());
        std::uniform_int_distribution<int64_t> local_dist(0, P - 1);
        int64_t local_vars[N_VARS] = {};
        int64_t local_t = 0;

        while (found.load(std::memory_order_relaxed) < target) {
            trials.fetch_add(1, std::memory_order_relaxed);

            if (slice_mode) {
                local_t = local_dist(local_rng);
                for (int k = 0; k < N_VARS; k++)
                    local_vars[k] = addmod(vec_a[k], mulmod(local_t, vec_b[k]));
            } else {
                for (int i = 0; i < N_VARS; i++)
                    local_vars[i] = local_dist(local_rng);
            }

            bool ok = true;
            for (const Node* r : check_order) {
                int64_t v;
                if (!eval_node(r, local_vars, v) || !is_square(v)) { ok = false; break; }
            }
            if (!ok) continue;

            // Claim a slot; discard if we've already hit the target
            if (found.fetch_add(1, std::memory_order_relaxed) >= target) break;

            std::string line;
            for (int i = 0; i < N_VARS; i++) {
                if (i) line += ' ';
                line += std::to_string(local_vars[i]);
            }
            if (slice_mode) line += ' ' + std::to_string(local_t);
            if (do_sqrts)
                for (const auto& r : roots) {
                    int64_t v = 0;
                    eval_node(r.get(), local_vars, v); // guaranteed to succeed: already checked above
                    line += ' ' + std::to_string(sqrtmod(v));
                }
            writeln(line);
        }
    }

    int64_t f = found.load();
    int64_t t = trials.load();
    if (f > target) f = target; // clamp display value
    std::cerr << "Found " << f << " in " << t << " trials"
              << " (" << (100.0 * f / t) << "% hit rate)\n";
    return 0;
}
