"""Write kaviar variants to only cgi, only illumina, and both"""
import csv
import sys

def writeInfo(dataSource):
    """A function to write information to a vcf file."""

    with open(dataSource, 'a') as ofh:
        rowWrite = '\t'.join(row) + '\n'
        ofh.write(rowWrite)


#Establish files to write variants to after classification/sorting
infile, outdir = sys.argv[1], sys.argv[2]
gsOutfile = outdir + 'good_sources.vcf'
cgiOutfile = outdir + 'cgi_only_sources.vcf'
illuminaOutfile = outdir + 'illumina_only_sources.vcf'
cgiIlluminaOutfile = outdir + 'cgi_illumina_sources.vcf'


with open(infile, 'r') as fh:  #open the vcf file in read mode specified

    readFile = csv.reader(fh, delimiter='\t')  #read the vcf file using csv module and specifying tab delimitation

    for row in readFile:

        if row[0][0] == '#':    #write heading lines of vcf read file to each output vcf file
            writeInfo(cgiIlluminaOutfile)
            writeInfo(cgiOutfile)
            writeInfo(illuminaOutfile)
            writeInfo(gsOutfile)

        else:
            info = row[7].split(';')  #split the info section of vcf file by semicolon to read data source

            for x in info:
                if x[0:2] == 'DS':    #check each entry in info section to find data source information
                    sources = info[info.index(x)].split('|')    #split the data source info by pipe and sort based on following conditions

                    goodSources = [x for x in sources if not x[0:2] in ('NA', 'GS') and not 'Wellderly' in x and x not in ('ISB_founders-Nge3', 'Inova_CGI_founders-Nge3')]

                    cgi_ls = [x for x in sources if x[0:2] in ('NA','GS') or 'Wellderly' in x or x in ('ISB_founders-Nge3','Inova_CGI_founders-Nge3')]

                    ill_ls = [x for x in sources if not x[0:2] in ('NA','GS') and not 'Wellderly' in x and not x in ('ISB_founders-Nge3','Inova_Illumina_founders-Nge3')]

                    #write the read data to corresponding vcf file calling the writeInfo function

                    if cgi_ls and ill_ls:
                        writeInfo(cgiIlluminaOutfile)


                    elif cgi_ls:
                        writeInfo(cgiOutfile)


                    elif ill_ls:
                        writeInfo(illuminaOutfile)


                    else:
                        writeInfo(gsOutfile)


                else:
                    pass

