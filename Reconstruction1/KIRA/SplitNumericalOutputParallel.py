
import sys, getopt, os
import multiprocessing as mp


def parse_args(argv):
    inputfile = ''
    outputfile = ''
    njobs = os.cpu_count()
    try:
        opts, args = getopt.getopt(argv, "hi:o:j:", ["ifile=", "ofile=", "jobs="])
    except getopt.GetoptError:
        print('SplitNumericalOutputParallel.py -i <inputfile> -o <outputfile> [-j <njobs>]')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('SplitNumericalOutputParallel.py -i <inputfile> -o <outputfile> [-j <njobs>]')
            sys.exit()
        if opt in ("-i", "--ifile"):
            inputfile = arg
        if opt in ("-o", "--ofile"):
            outputfile = arg
        if opt in ("-j", "--jobs"):
            njobs = int(arg)
    return inputfile, outputfile, njobs


def split_one_file(task):
    filePos, inputfile, outputfile = task
    inputfileNew = inputfile[:-3] + str(filePos) + ".m"

    print(f"[filePos {filePos}] Input file is {1000*filePos+1}")
    outFolder = f"{filePos:04d}"
    os.makedirs(f"{outputfile}/{outFolder}", exist_ok=True)
    linecount = 1000 * filePos + 1
    file1 = open(f"{outputfile}/{outFolder}/phasespace_{linecount:07d}.m", 'w')
    fileR = open(inputfileNew, 'r')
    prevline = ""

    linecount = linecount - 1

    for i, line in enumerate(fileR):
        if i == 0:
            primeline = line
            file1.writelines(primeline)

        if (not "{HTP1" in line) and ("{{" in line):
            linecount += 1
            if linecount > 1:
                file1.writelines(prevline[:-2] + "}")
                file1.close()
                file1 = open(f'{outputfile}/{outFolder}/phasespace_{linecount:07d}.m', 'w')
                print(f"[filePos {filePos}] Current SamplePoint: {linecount:07d}")
                file1.writelines(primeline)
        elif i > 0:
            if "}}}}," in prevline:
                prevline = prevline[:-2]
            file1.writelines(prevline)
        prevline = line

    file1.writelines(prevline)
    file1.close()
    fileR.close()
    return filePos


if __name__ == "__main__":
    inputfile, outputfile, njobs = parse_args(sys.argv[1:])
    os.makedirs(outputfile, exist_ok=True)

    # Each raw Kira output file (filePos 0..999) is split independently into
    # its own output subfolder, so different filePos values can be processed
    # in parallel with no shared state between them.
    tasks = []
    for filePos in range(1000):
        inputfileNew = inputfile[:-3] + str(filePos) + ".m"
        if os.path.isfile(inputfileNew):
            tasks.append((filePos, inputfile, outputfile))

    print(f"Found {len(tasks)} input file(s) to split, using {njobs} worker process(es)")

    with mp.Pool(processes=njobs) as pool:
        for result in pool.imap_unordered(split_one_file, tasks):
            print(f"Finished filePos {result}")
