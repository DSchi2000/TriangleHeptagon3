
import sys, getopt, os
def main(argv):
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print ('SplitNumericalOutput.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('SplitNumericalOutput.py -i <inputfile> -o <outputfile>')
            sys.exit()
        if opt in ("-i", "--ifile"):
            inputfile = arg
        if opt in ("-o", "--ofile"):
            outputfile = arg
    return (inputfile,outputfile)

inputfile = ''
outputfile = ''
if __name__ == "__main__":
   inputfile,outputfile=main(sys.argv[1:])

os.system(f"mkdir -p {outputfile}/")
for filePos in range(1000):
    inputfileNew=inputfile[:-3]+str(filePos)+".m"
    if(os.path.isfile(inputfileNew)):
        print ('Input file is ', 1000*filePos+1)
        outFolder=f"{filePos:04d}"
        os.system(f"mkdir -p {outputfile}/{outFolder}")
        linecount=1000*filePos+1
        file1 = open(f"{outputfile}/{outFolder}/phasespace_{linecount:07d}.m", 'w')
        fileR = open(inputfileNew, 'r') 
        prevline=""

        linecount=linecount-1

        for i,line in enumerate(fileR):
            if i==0:
                primeline=line
                file1.writelines(primeline)

            if (not "{HTP1" in line) and ("{{" in line):
                linecount+=1
                # print(linecount)
                if (linecount>1):
                    file1.writelines(prevline[:-2]+"}")
                    file1.close()
                    file1 = open(f'{outputfile}/{outFolder}/phasespace_{linecount:07d}.m', 'w')
                    print(f"Current SamplePoint: {linecount:07d}")
                    file1.writelines(primeline)
            elif i>0:
                if "}}}}," in prevline:
                    prevline=prevline[:-2]
                file1.writelines(prevline)
            prevline=line

        file1.writelines(prevline)
        file1.close()










