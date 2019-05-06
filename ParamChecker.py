#!/usr/bin/env python
import logging
import os


class ParamChecker:
    # Constructor that creates two empty arrays that
    def __init__(self):
        # ???: Maybe rename this __name__ to track the log outputs easier?
        self.vaseLogger = logging.getLogger("VaSe_Logger")
        self.vcfFolders = []
        self.bamFolders = []
        self.acceptorBam = ""
        self.fastqIn1 = ""
        self.fastqIn2 = ""
        self.outDir = ""
        self.fastqOutLocation = ""
        self.varConOutLocation = ""
        self.donorBreadOutLocation = ""
        self.acceptorBreadOutLocation = ""
        self.logLocation = ""

    # Check the logging parameter to determine where to
    # write the logfile to.
    # ???: Can we grab the current UUID of the VaSe run to assign a
    # default, unique log name? e.g. VaSeBuilder.o123456.log or something
    # Also, I think this can be defined as an 'action' in argparse, so
    # that it is implemented directly in the parsed argument list.
    # Also, it doesn't seem like there is currently a way to handle when
    # a user specifies a --log argument that is not a current file, is
    # not a directory, but also does not end in 'log' or 'txt'.

# =============================================================================
#     def checkLog(self, log_arg):

#         if (log_arg is None):
#             self.logLocation = 'VaSeBuilder.log'  # or some unique default? ^^^
#         elif (os.path.isdir(log_arg)):
#             self.logLocation = log_arg + "/VaSeBuilder.log"
#         elif (not (os.path.isfile(log_arg))
#               and log_arg.endswith((".log", ".txt."))):
#             self.logLocation = log_arg
#         elif (not (os.path.isfile(log_arg))
#               and not log_arg.endswith((".log", ".txt."))):
#             self.logLocation = log_arg + ".log"

#         return self.logLocation
# =============================================================================

    def checkLog(self, logParam):
        logloc = "VaSeBuilder.log"

        if (logParam is not None):
            # Check the location of the log file if the --log parameter
            # has been set.
            if (not(os.path.isfile(logParam))
               and (logParam.endswith(".log") or logParam.endswith(".txt"))):
                logloc = logParam

            # Check to make sure the provided --log parameter value is
            # not a directory (Directories could be named
            # 'something.log')
            if (os.path.isdir(logParam)):
                logloc = logParam + "/VaSeBuilder.log"
        self.logLocation = logloc
        return self.logLocation

    # Checks whether provided folders exist
    # ???: I'm not really clear on how this is checking if folders exist.
    # I feel like there might be a simpler way to check this. 
    # Maybe use os.walk to skip making a
    # list of folders to instead make a list of available VCF files in
    # the directories provided.
    def checkFoldersExist(self, paramVals, fileExt):
        existingFolders = []

        # Check whether the provided folders exist.
        for parval in paramVals:
            foldername = self.getFolderName(parval)
            if (foldername != parval):
                self.vaseLogger.warning("File instead of folder was supplied. "
                                        + "Using the folder "
                                        + foldername
                                        + " as input folder")

            # Check if the supplied value is a folder or not and
            # contains any vcf/bam files.
            if (not (os.path.isdir(foldername))):
                self.vaseLogger.warning("Folder "
                                        + foldername +
                                        + " was not found and will"
                                        + " therefore be skipped")
            else:
                if (self.checkFolderContents(foldername, fileExt) == 0):
                    self.vaseLogger.warning("Folder "
                                            + foldername
                                            + " exists but contains no "
                                            + fileExt
                                            + " files")
                else:
                    self.vaseLogger.info("Folder "
                                         + foldername
                                         + " will be included")
                    existingFolders.append(foldername)
        return existingFolders

    # Checks whether at least one file with a provided extension
    # (.vcf or .bam) is present.
    def checkFolderContents(self, folderToCheck, fileExt):
        vbCnt = 0
        for vbFile in os.listdir(folderToCheck):
            if (vbFile.endswith("." + fileExt)):
                vbCnt += 1
        self.vaseLogger.debug("Folder "
                              + folderToCheck
                              + " contains "
                              + str(vbCnt)
                              + " "
                              + fileExt
                              + " files")
        return vbCnt

    # Checks whether a provided file exists.
    def checkFileExists(self, fileLoc):
        if (os.path.isfile(fileLoc)):
            self.vaseLogger.debug("File " + fileLoc + " exists")
            return True
        self.vaseLogger.debug("File " + fileLoc + " does not exist")
        return False

    # Return the directory name of an output location.
    def isValidOutputLocation(self, outFileName):
        return (os.path.isdir(os.path.dirname(outFileName)))

    # Checks whether the values of the parameters are correct
    # (do files/folders exist for example)
    # [Function should perhaps be split into smaller functions]
    # XXX: I think a lot of this is built-in with argparse. That would
    # probably already avoid having to loop over the arglist and having
    # to check the identity of each one to apply the proper test. If need
    # be, we can probably split these tests for each argument and pass
    # them to argparse (I think)
    def checkParameters(self, vaseArgVals):

        # Loop over the provided parameters.
        for param in vaseArgVals:

            # If the current parameter is vcfin, check whether there are
            # any valid VCF folders to use.
            if (param == 'donorvcf'):
                vcfFolders = self.checkFoldersExist(vaseArgVals[param],
                                                    "vcf.gz")
                # XXX: This is probably unnecessary if '+' is specified in addarg
                if (len(vcfFolders) == 0):
                    self.vaseLogger.critical(
                            "No folders containing VCF files were found. "
                            "Please supply existing folders next time :)"
                            )
                    return False
                self.vcfFolders = vcfFolders

            # If the current parameter is bamin, check whether there are
            # any valid BAM folders to use.
            if (param == 'donorbam'):
                bamFolders = self.checkFoldersExist(vaseArgVals[param], "bam")
                # XXX: This is probably unnecessary if '+' is specified in addarg
                if (len(bamFolders) == 0):
                    self.vaseLogger.critical(
                            "No folders containing BAM files were found. "
                            "Please supply existing folders next time :)"
                            )
                    return False
                self.bamFolders = bamFolders

            # If the current parameter is bam, check whether a valid BAM
            # file is provided.
            if (param == 'acceptorbam'):
                if (not self.checkFileExists(vaseArgVals[param])):
                    self.vaseLogger.critical("No valid NIST BAM file "
                                             "supplied :(")
                    return False
                self.acceptorBam = vaseArgVals[param]

            # If the current parameter is valfastq1, check whether a
            # valid R1 fastq file is provided.
            if (param == 'templatefq1'):
                if (not self.checkFileExists(vaseArgVals[param])):
                    self.vaseLogger.critical("Provided R1 FastQ input file "
                                             "does not exist")
                    return False
                self.fastqIn1 = vaseArgVals[param]

            # If the current parameter is valfastq2, check whether a
            # valid R2 fastq file is provided.
            if (param == 'templatefq2'):
                if (not self.checkFileExists(vaseArgVals[param])):
                    self.vaseLogger.critical("Provided R2 FastQ input file "
                                             "does not exist")
                    return False
                self.fastqIn2 = vaseArgVals[param]

            # If the current parameter is out, check whether it is a
            # valid output location.
            if (param == 'out'):
                if (not self.isValidOutputLocation(vaseArgVals[param])):
                    return False
                self.outDir = vaseArgVals[param]

            # If the current parameters is fastqout, check if a name has
            # been provided
            if (param == 'fastqout'):
                self.fastqOutLocation = self.getOutputName(vaseArgVals[param],
                                                           "VaSe")

            # If the current parameter is varcon, check whether a valid
            # output location is provided
            if (param == 'varcon'):
                self.varConOutLocation = self.getOutputName(vaseArgVals[param],
                                                            'varcon.txt')

            # If the current parameters is varbread, check whether a
            # valid output location is provided
            if (param == 'donorbread'):
                self.donorBreadOutLocation = self.getOutputName(
                        vaseArgVals[param],
                        'donorbread.txt'
                        )

            # If the current parameter is nistbread, check whether a
            # valid output location is provided
            if (param == 'acceptorbread'):
                self.acceptorBreadOutLocation = self.getOutputName(
                        vaseArgVals[param],
                        'acceptorbread.txt'
                        )

        # Return the lists of valid VCF and BAM folders that can be used
        # by the program.
        return True

    # Returns thename of the folder name of a parameter value (if the
    # parameter value is )
    # XXX: I don't really think this needs to be its own fxn
    def getFolderName(self, foldername):
        if (os.path.isfile(foldername) or (not os.path.isdir(foldername))):
            return os.path.dirname(foldername)
        return foldername

    # Returns the name of an output file (is used for parameters
    # fastqout, varcon, donorbread and acceptorbread)
    def getOutputName(self, outfilename, defaultoutname):
        if (outfilename is not None):
            if ("/" in outfilename):
                return outfilename.split("/")[-1]
            elif ("\\" in outfilename):
                return outfilename.split("\\")[-1]
            return outfilename
        return defaultoutname

    # Returns the list of valid VCF folders.
    def getValidVcfFolders(self):
        return self.vcfFolders

    # Returns the list of valid BAM folders.
    def getValidBamFolders(self):
        return self.bamFolders

    # Returns the location of the  NIST BAM file.
    def getAcceptorBam(self):
        return self.acceptorBam

    # Returns the location and name of the first (R1) fastq input file.
    def getFirstFastqInLocation(self):
        return self.fastqIn1

    # Returns the location and name of the second (R2) fastq input file.
    def getSecondFastqInLocation(self):
        return self.fastqIn2

    # Returns the location(s) and names of the two (R1 and R2) fastq
    # input files.
    def getFastqInLocations(self):
        return [self.fastqIn1, self.fastqIn2]

    # Returns the location to write the output to.
    def getOutDirLocation(self):
        return self.outDir

    # Returns the location of the FastQ file that will be produced by
    # VaSeBuilder.
    def getFastqOutLocation(self):
        return self.outDir+"/"+self.fastqOutLocation

    # Returns the location of file that will contain the variants and
    # their context start and stops.
    def getVariantContextOutLocation(self):
        return self.outDir+"/"+self.varConOutLocation

    # Returns the location of the file that will contain the variants
    # and their associatied patient BAM reads.
    def getDonorBamReadOutLocation(self):
        return self.outDir+"/"+self.donorBreadOutLocation

    # Returns the location of the file that will containt the variants
    # and their associated NIST reads.
    def getAcceptorBamReadOutLocation(self):
        return self.outDir+"/"+self.acceptorBreadOutLocation

    # Retuns the location to write the log file(s) to.
    def getLogFileLocation(self):
        return self.logLocation
