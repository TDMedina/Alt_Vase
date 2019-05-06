#!/usr/bin/env python
import logging
import os
import pysam


class VcfBamScanner:

    # Constructor that creates two empty hashmaps (dicts)
    def __init__(self):
        # XXX: Maybe change the logger name to '__name__' so we can
        # track the log easier?
        self.vaseLogger = logging.getLogger("VaSe_Logger")
        self.vcfSampleMap = {}
        self.bamSampleMap = {}

    # Scans the folders containing VCF files and returns a map that
    # links sample ids with vcf files.
    # XXX: Maybe we skip this whole scanning fxn by using some kind of sys find
    # thing instead?
    # XXX: Can we combine the VCF/BAM scanning functions? They seem to
    # perform nearly identical functions. For their differences,
    # a few IFs can be thrown in that make the decision based on an
    # additional keyword arg in the definition call.
    # XXX: def scanFolders(self, file_folders, file_ext):
    def scanVcfFolders(self, vcfFolders):
        # XXX: self.vaseLogger.info(f"Start scanning {file_ext.upper()} files.")
        # Do this for other log strings.
        self.vaseLogger.info("Start scanning VCF files")

        for vcfFolder in vcfFolders:
            # ???: Can we un-nest the IF statement below with a
            # guard clause instead? Might be tidier.
            if (os.path.isdir(vcfFolder)):
                self.vaseLogger.info("Scanning VCF files in " + vcfFolder)
                # XXX: Maybe include the IF clause in the FOR statement
                # to reduce nesting? Generator statement?
                # gen = (vcf for vcf in os.listdir(vcffolder)
                #        if vcf.endswith(('.vcf','.vcf.gz')))
                # for vcf in gen:
                for vcfFileName in os.listdir(vcfFolder):
                    # XXX: Maybe another guard clause here?
                    # XXX: if not (vcfFileName.endswith(f'.{file_ext}', f'.{file_ext}.gz'))):
                    #          pass
                    # XXX: if (vcfFileName.endswith(f'.{file_ext}', f'.{file_ext}.gz'))):
                    if (vcfFileName.endswith(".vcf")
                       or vcfFileName.endswith(".vcf.gz")):
                        # Do something to scan the file for the sample
                        try:
                            vcfFile = pysam.VariantFile(
                                    vcfFolder+"/"+vcfFileName,
                                    'r'
                                    )
                            self.vaseLogger.debug("Scanning VCF file "
                                                  + vcfFolder
                                                  + "/"
                                                  + vcfFileName)
                            # ???: Is there only one sample per VCF?
                            if (len(vcfFile.header.samples) > 0):
                                # This is the sample identifier
                                sampleid = vcfFile.header.samples[0]
                                self.vcfSampleMap[sampleid] = (vcfFolder
                                                               + "/"
                                                               + vcfFileName)
                            vcfFile.close()

                        except IOError as ioe:
                            self.vaseLogger.warning("VCF file "
                                                    + vcfFolder
                                                    + "/"
                                                    + vcfFileName
                                                    + " does not exist")
                        # ???: Is the IOError likely to occur? The file
                        # list is taken directly from OS, so would they
                        # have moved? Do other exceptions need to be
                        # taken into account?
            else:
                self.vaseLogger.info("Folder "
                                     + vcfFolder
                                     + " does not exist.")
        self.vaseLogger.info("Finished scanning VCF files")
        return self.vcfSampleMap

    # Scans the folders containing BAM files and returns a map that links
    # sample ids with bam files.
    def scanBamFolders(self, bamFolders):
        self.vaseLogger.info("Start scanning BAM files")

        # Scan BAM files in all provided folders
        for bamFolder in bamFolders:
            # XXX: Make this a guard clause maybe.
            if (os.path.isdir(bamFolder)):
                self.vaseLogger.info("Scanning BAM files in " + bamFolder)
                # XXX: Maybe use a generator here:
                # bam_gen = (bam for bam in os.listdir(bamFolder)
                #            if bam.endswith('bam'))
                # for bam in bam_gen:
                for bamFileName in os.listdir(bamFolder):
                    if (bamFileName.endswith(".bam")):
                        try:
                            bamFile = pysam.AlignmentFile(
                                    bamFolder+"/"+bamFileName,
                                    'rb'
                                    )
                            self.vaseLogger.debug("Scanning BAM file "
                                                  + bamFolder
                                                  + "/"
                                                  + bamFileName)
                            if (self.bamHasSampleName(bamFile)):
                                # The sample identifier
                                sampleid = bamFile.header["RG"][0]["SM"]
                                self.bamSampleMap[sampleid] = (bamFolder
                                                               + "/"
                                                               + bamFileName)

                            bamFile.close()

                        except IOError as ioe:
                            self.vaseLogger.warning("BAM file "
                                                    + bamFolder
                                                    + "/"
                                                    + bamFileName
                                                    + " does not exist")
            else:
                self.vaseLogger.info("Folder " + bamFolder + "does not exist.")

        self.vaseLogger.info("Finished scanning BAM files")
        return self.bamSampleMap

    # Checks whether the BAM file contains a sample name.
    def bamHasSampleName(self, bamFile):
        if ("RG" in bamFile.header):
            if (len(bamFile.header["RG"]) > 0):
                if ("SM" in bamFile.header["RG"][0]):
                    return True
        return False

# =============================================================================
# ^^^^def bamHasSampleName(self, bamFile):
#         if ('RG' not in bamFile.header):
#             return False
#         if (len(bamFile.header['RG']) == 0):
#             return False
#         if ('SM' not in bamFile.header['RG'][0]):
#             return False
#         else:
#             return True
# =============================================================================

    # Returns the map that links VCF files and samples
    def getVcfSampleMap(self):
        return self.vcfSampleMap

    # Returns the map that links BAM files and samples
    def getBamSampleMap(self):
        return self.bamSampleMap

    # Returns a map that links VCF files to BAM files
    # ???: I'm not sure that this actually does anything. It checks if
    # any ID's are found in both the VCF and BAM map dicts, then makes
    # a dict of the matching IDs' paths mapping to each other. I think
    # this is an extra unecessary list that doesn't add any new info, and also
    # has very long keys. 
    #
    # If the point is to make
    # sure that each VCF has a matching BAM file, I think a list of
    # ID's (keys) found in both dictionaries would make more sense.
    # {'A': 'john.vcf', 'B': 'jack.vcf', 'C': 'jill.vcf'}
    # {'A': 'john.bam', 'C': 'jill.bam'}
    # ['A', 'C']
    # Also would this be better as part of the class definition, so that the
    # VCFtoBam map (list of common keys) is part of the class, and the method
    # simply returns it? This might be more consistent.
    # List version:
# =============================================================================
#     def getVcfToBamMap(self):
#         vcfToBamMap = [vcf_id for vcf_id in self.vcfSampleMap
#                        if vcf_id in self.bamSampleMap]
#         return vcfToBamMap
# =============================================================================
    def getVcfToBamMap(self):
        vcfToBamMap = {}
        for sampleid in self.vcfSampleMap:
            if (sampleid in self.bamSampleMap):
                vcfToBamMap[self.vcfSampleMap[sampleid]] = self.bamSampleMap[sampleid]
        return vcfToBamMap
