# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 17:12:20 2019

@author: medinatd
"""

from alt_folder_scan_fxn import file_mapper


def main(self):
    vaseArgList = self.getVaSeParameters()
    pmc = ParamChecker()
    self.vaseLogger = self.startLogger(pmc, vaseArgList['log'])

    if (not (pmc.checkParameters(vaseArgList))):
        self.vaseLogger.critical("Not all parameters are correct. "
                                 "Please check log for more info.")
        exit()

    file_maps = file_mapper()
    file_maps.file_scan('bam', root_folders=pmc.getValidBamFolders())
    file_maps.file_scan('vcf', root_folders=pmc.getValidVcfFolders())
    file_maps.get_common_ids()

    vaseB = VaSeBuilder(uuid.uuid4().hex)

    if (len(file_maps.common_ids == 0)):
        self.vaseLogger.critical("No valid samples available to "
                                 "create new validation set")
        exit()

    # Start the procedure to build the validation set.
    vaseB.buildValidationSet(file_maps.common_ids,
                             file_maps.vcf_map,
                             file_maps.bam_map,
                             pmc.getAcceptorBam(),
                             pmc.getFirstFastqInLocation(),
                             pmc.getSecondFastqInLocation(),
                             pmc.getOutDirLocation(),
                             pmc.getFastqOutLocation(),
                             pmc.getVariantContextOutLocation(),
                             pmc.getDonorBamReadOutLocation(),
                             pmc.getAcceptorBamReadOutLocation())
    self.vaseLogger.info("VaSeBuilder run completed successfully.")
