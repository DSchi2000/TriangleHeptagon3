// 6-point phase space generator over GF(p), p = 2^31-1
//
// Reads Roots.txt  — one polynomial per line; all must be QR mod p.
// Reads SixPointGram.txt — quadratic constraint that determines the aux variable.
//
// Algorithm:
//   1. Sample free variables uniformly in GF(p).
//   2. Extract A, B, C s.t. Gram = A*aux^2 + B*aux + C, solve for aux.
//   3. For each aux solution, check all root conditions via Jacobi symbol.
//   4. Output valid points.
//
// Adding/removing roots: just edit Roots.txt — no recompilation needed.

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
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
    "d","mt2","s12","s123","s16","s23","s234","s34","s345","s45","s56"
};
static constexpr int N_VARS = 11;

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

// ── Polynomial representation ─────────────────────────────────────────────────
struct Monomial {
    int64_t coeff;
    int8_t  exp[N_VARS];
};
using Poly = std::vector<Monomial>;

static int64_t eval_poly(const Poly& poly, const int64_t* vars) {
    int64_t sum = 0;
    for (const auto& m : poly) {
        int64_t t = m.coeff;
        for (int i = 0; i < N_VARS; i++)
            if (m.exp[i]) t = mulmod(t, powmod(vars[i], m.exp[i]));
        sum = addmod(sum, t);
    }
    return sum;
}

// ── Parser ────────────────────────────────────────────────────────────────────
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
    char get()   { skip(); return p < end ? *p++ : 0; }
    bool eat(char c) { if (peek()==c) { get(); return true; } return false; }
    bool at_end()   { skip(); return p >= end; }

    int64_t read_int() {
        skip(); int64_t v = 0;
        while (p < end && isdigit((unsigned char)*p)) v = v * 10 + (*p++ - '0');
        return v;
    }
    std::pair<const char*, int> read_ident() {
        skip(); const char* s = p;
        while (p < end && (isalnum((unsigned char)*p) || *p == '_')) ++p;
        return {s, (int)(p - s)};
    }
};

static Monomial parse_monomial(Lex& lx) {
    Monomial m; m.coeff = 1; memset(m.exp, 0, sizeof(m.exp));

    bool in_paren = lx.eat('(');
    bool first = true;

    while (true) {
        char c = lx.peek();
        if (!in_paren && (c == 0 || c == '+' || c == '-')) break;
        if (c == ')') { if (in_paren) { lx.get(); in_paren = false; } break; }
        if (!first && !lx.eat('*')) break;
        first = false;

        c = lx.peek();
        if (isdigit((unsigned char)c)) {
            m.coeff = mulmod(m.coeff, lx.read_int() % P);
        } else if (isalpha((unsigned char)c)) {
            auto [s, len] = lx.read_ident();
            int idx = var_index(s, len);
            int e = 1;
            if (lx.eat('^')) e = (int)lx.read_int();
            m.exp[idx] += (int8_t)e;
        } else break;
    }

    if (lx.eat('/')) {
        int64_t d = lx.read_int();
        m.coeff = mulmod(m.coeff, invmod(d % P));
    }
    return m;
}

static Poly parse_poly(const char* s, size_t n) {
    Poly result;
    Lex lx{s, s + n};
    while (!lx.at_end()) {
        int sign = 1;
        char c = lx.peek();
        if      (c == '+') { lx.get(); sign =  1; }
        else if (c == '-') { lx.get(); sign = -1; }
        if (lx.at_end()) break;

        Monomial m = parse_monomial(lx);
        if (sign < 0) m.coeff = submod(0, m.coeff);
        if (m.coeff) result.push_back(m);
    }
    return result;
}

static Poly parse_line(const std::string& line) {
    const char* p = line.c_str();
    size_t n = line.size();
    const char* q = p;
    while (q - p < 12 && q < p + n && isdigit((unsigned char)*q)) ++q;
    if (q > p && q < p + n && *q == '\t') { p = q + 1; n = line.size() - (p - line.c_str()); }
    return parse_poly(p, n);
}

static std::vector<Poly> load_roots(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("Cannot open: " + path);
    std::vector<Poly> roots;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        auto p = parse_line(line);
        if (!p.empty()) roots.push_back(std::move(p));
    }
    return roots;
}

static Poly load_gram(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("Cannot open: " + path);
    std::string content((std::istreambuf_iterator<char>(f)),
                         std::istreambuf_iterator<char>());
    return parse_line(content);
}

// ── Quadratic solver mod P ────────────────────────────────────────────────────
static int solve_quad(int64_t A, int64_t B, int64_t C,
                       int64_t& x1, int64_t& x2) {
    if (A == 0) {
        if (B == 0) { return C == 0 ? -1 : 0; }
        x1 = mulmod(submod(0, C), invmod(B));
        return 1;
    }
    int64_t D = submod(mulmod(B, B), mulmod(4, mulmod(A, C)));
    int jD = jacobi(D, P);
    if (jD < 0) return 0;
    int64_t sqD   = sqrtmod(D);
    int64_t inv2A = invmod(mulmod(2, A));
    x1 = mulmod(submod(submod(0, B), sqD), inv2A);
    x2 = mulmod(addmod(submod(0, B), sqD), inv2A);
    return jD == 0 ? 1 : 2;
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
        << "  --roots   <path>    Roots.txt          (default: Roots.txt)\n"
        << "  --gram    <path>    SixPointGram.txt   (default: SixPointGram.txt)\n"
        << "  --auxvar  <name>    variable solved from Gram (default: s16)\n"
        << "  --npoints <N>       points to generate  (default: 100)\n"
        << "  --seed    <S>       RNG seed             (default: 42)\n"
        << "  --threads <N>       number of CPU cores  (default: 1)\n"
        << "  --out     <file>    also save to file\n"
        << "  --sqrt              also output square roots of each expression\n"
        << "  --slice <a> <b>     slice mode: s = a + t*b for random t\n"
        << "  --gram-resolved     aux variable is already fixed by the slice\n"
        << "                      (e.g. a BCFW slice already satisfies G6=0) —\n"
        << "                      don't solve the Gram quadratic again; requires\n"
        << "                      --slice with vectors covering all " << N_VARS << " variables\n"
        << "  N                   positional: equivalent to --npoints N\n";
}

// ── Main ─────────────────────────────────────────────────────────────────────
int main(int argc, char* argv[]) {
    std::string roots_path   = "Roots.txt";
    std::string gram_path    = "SixPointGram.txt";
    std::string auxvar_name  = "s16";
    int64_t     target       = 100;
    uint64_t    seed         = std::random_device{}();
    int         n_threads    = 1;
    bool        do_sqrts     = false;
    bool        gram_resolved = false;
    std::string out_path;
    std::string slice_a_path, slice_b_path;

    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        if      (a == "--roots"   && i+1 < argc) roots_path   = argv[++i];
        else if (a == "--gram"    && i+1 < argc) gram_path    = argv[++i];
        else if (a == "--auxvar"  && i+1 < argc) auxvar_name  = argv[++i];
        else if (a == "--out"     && i+1 < argc) out_path     = argv[++i];
        else if (a == "--npoints" && i+1 < argc) target       = std::stoll(argv[++i]);
        else if (a == "--seed"    && i+1 < argc) seed         = std::stoull(argv[++i]);
        else if (a == "--threads" && i+1 < argc) n_threads    = std::stoi(argv[++i]);
        else if (a == "--slice"   && i+2 < argc) { slice_a_path = argv[++i]; slice_b_path = argv[++i]; }
        else if (a == "--sqrt")  do_sqrts = true;
        else if (a == "--gram-resolved") gram_resolved = true;
        else if (a == "--help")  { usage(argv[0]); return 0; }
        else if (isdigit(a[0])) target = std::stoll(a);
    }

    // Resolve aux variable index
    int idx_aux = -1;
    for (int i = 0; i < N_VARS; i++)
        if (auxvar_name == VAR_NAMES[i]) { idx_aux = i; break; }
    if (idx_aux < 0) {
        std::cerr << "Error: unknown --auxvar '" << auxvar_name << "'\n";
        usage(argv[0]); return 1;
    }

    // Build free_vars = all indices except idx_aux
    std::vector<int> free_vars;
    for (int i = 0; i < N_VARS; i++)
        if (i != idx_aux) free_vars.push_back(i);

    bool slice_mode = !slice_a_path.empty();
    if (gram_resolved && !slice_mode) {
        std::cerr << "Error: --gram-resolved requires --slice "
                      "(the aux variable must come from the slice)\n";
        return 1;
    }

    // slice_idx lists which variables are driven by the slice line s = a + t*b.
    // Normally that's the free variables only, and the aux variable is solved
    // from the Gram constraint. With --gram-resolved the slice already fixes
    // the aux variable too (e.g. a BCFW slice that satisfies G6=0 identically),
    // so it is included here and the Gram quadratic is never solved.
    std::vector<int> slice_idx = free_vars;
    if (gram_resolved) {
        slice_idx.clear();
        for (int i = 0; i < N_VARS; i++) slice_idx.push_back(i);
    }
    int n_slice = (int)slice_idx.size();

    // Load slice vectors if needed
    std::vector<int64_t> vec_a, vec_b;
    if (slice_mode) {
        try {
            vec_a = load_vec(slice_a_path, n_slice);
            vec_b = load_vec(slice_b_path, n_slice);
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << "\n"; return 1;
        }
        std::cerr << "Slice mode: s = a + t*b  (t random in GF(p))\n";
        if (gram_resolved)
            std::cerr << "Gram constraint: already resolved by the slice "
                          "(not solved again)\n";
    }

    std::vector<Poly> roots;
    Poly gram;
    try {
        roots = load_roots(roots_path);
        if (!gram_resolved) gram = load_gram(gram_path);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n"; return 1;
    }
    std::cerr << "Seed: " << seed << "\n";
    std::cerr << "Aux variable: " << auxvar_name << "\n";
    std::cerr << "Loaded " << roots.size() << " root conditions";
    if (!gram_resolved) std::cerr << ", " << gram.size() << " gram monomials";
    std::cerr << "\n";

    // Split gram into gramC + aux*gramB + aux^2*gramA
    Poly gramA, gramB, gramC;
    if (!gram_resolved) {
        for (auto m : gram) {
            int e = m.exp[idx_aux];
            m.exp[idx_aux] = 0;
            if      (e == 0) gramC.push_back(m);
            else if (e == 1) gramB.push_back(m);
            else if (e == 2) gramA.push_back(m);
            else {
                std::cerr << "Error: Gram has " << auxvar_name << "^" << e
                          << " > 2 — choose a different aux variable\n";
                return 1;
            }
        }
    }

    // Split roots into pre-gram (no aux) and post-gram (has aux)
    std::vector<Poly> pre_roots, post_roots;
    for (auto& r : roots) {
        bool has_aux = false;
        for (const auto& m : r) if (m.exp[idx_aux]) { has_aux = true; break; }
        (has_aux ? post_roots : pre_roots).push_back(r);
    }
    auto by_size = [](const Poly& a, const Poly& b){ return a.size() < b.size(); };
    std::sort(pre_roots.begin(),  pre_roots.end(),  by_size);
    std::sort(post_roots.begin(), post_roots.end(), by_size);
    std::cerr << "  " << pre_roots.size()  << " pre-gram roots, "
              << post_roots.size() << " post-gram roots\n";
    std::cerr << "Threads: " << n_threads << "\n";

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
        // Each thread gets its own RNG seeded uniquely
        std::mt19937_64 local_rng(seed + (uint64_t)omp_get_thread_num());
        std::uniform_int_distribution<int64_t> local_dist(0, P - 1);
        int64_t local_vars[N_VARS] = {};
        int64_t local_t = 0;

        while (found.load(std::memory_order_relaxed) < target) {
            trials.fetch_add(1, std::memory_order_relaxed);

            // Sample free variables (slice_idx also covers the aux variable
            // when --gram-resolved is set)
            if (slice_mode) {
                local_t = local_dist(local_rng);
                for (int k = 0; k < n_slice; k++)
                    local_vars[slice_idx[k]] = addmod(vec_a[k], mulmod(local_t, vec_b[k]));
            } else {
                for (int i = 0; i < N_VARS; i++)
                    if (i != idx_aux) local_vars[i] = local_dist(local_rng);
            }

            // Pre-gram check
            bool pre_ok = true;
            for (const auto& r : pre_roots)
                if (!is_square(eval_poly(r, local_vars))) { pre_ok = false; break; }
            if (!pre_ok) continue;

            // Solve Gram for aux variable (skipped when --gram-resolved:
            // the slice already fixed local_vars[idx_aux] above)
            if (!gram_resolved) {
                int64_t C = eval_poly(gramC, local_vars);
                int64_t B = eval_poly(gramB, local_vars);
                int64_t A = eval_poly(gramA, local_vars);

                int64_t x1, x2;
                int nsol = solve_quad(A, B, C, x1, x2);
                if (nsol == 0) continue;

                if (nsol == -1) local_vars[idx_aux] = local_dist(local_rng);
                else            local_vars[idx_aux] = x1;
            }

            // Post-gram check
            bool post_ok = true;
            for (const auto& r : post_roots)
                if (!is_square(eval_poly(r, local_vars))) { post_ok = false; break; }
            if (!post_ok) continue;

            // Claim a slot; discard if we've already hit the target
            if (found.fetch_add(1, std::memory_order_relaxed) >= target) break;

            // Build and emit output line
            std::string line;
            for (int i = 0; i < N_VARS; i++) {
                if (i) line += ' ';
                line += std::to_string(local_vars[i]);
            }
            if (slice_mode) line += ' ' + std::to_string(local_t);
            if (do_sqrts)
                for (const auto& r : roots)
                    line += ' ' + std::to_string(sqrtmod(eval_poly(r, local_vars)));
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
