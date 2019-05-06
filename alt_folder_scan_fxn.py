# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 15:47:23 2019

@author: medinatd
"""

import logging
import os
import pysam


class file_mapper:

    # Constructor that creates two empty hashmaps (dicts)
    def __init__(self):
        self.vaseLogger = logging.getLogger('__name__')
        self.vcf_map = {}
        self.bam_map = {}
        self.common_ids = []

    def file_scan(self, file_ext, root_folders='.'):

        # XXX: Maybe add an optional flag to disable/enable recursive walk
        # through provided dirs
        self.vaseLogger.info(f'Searching for {file_ext.upper()} files...')
        found_files = []
        for root in list(root_folders):
            self.vaseLogger.debug(f'Searching in {root}...')
            if (not (os.path.isdir(root))):
                print(f'Directory "{root}" not found; Skipping.')
                continue
            tree = list(os.walk(root))
            found = [(y[0] + '\\' + x) for y in tree for x in y[2]
                     if x.endswith((f'.{file_ext}', f'.{file_ext}.gz'))]
            found_files = found_files + found

        found_files = list(set(found_files))
        if (file_ext == 'vcf'):
            for file_path in found:
                self.vaseLogger.debug(f'Scanning VCF file {file_path}...')
                try:
                    vcf_file = pysam.VariantFile(file_path)
                except IOError:
                    self.vaseLogger.warning('Unable to read file '
                                            f'{file_path}; Skipping file.')
                    continue
                except ModuleNotFoundError:
                    print('This is a test')
                if(len(vcf_file.header.samples) > 0):
                    vcf_id = vcf_file.header.samples[0]
                    self.vcf_map[vcf_id] = file_path
                vcf_file.close()
            self.vaseLogger.info('Finished scanning VCF files.')
            return self.vcf_map

        elif (file_ext == 'bam'):
            for file_path in found:
                self.vaseLogger.debug(f'Scanning BAM file {file_path}...')
                try:
                    bam_file = pysam.AlignmentFile(file_path)
                except IOError:
                    self.vaseLogger.warning('Unable to read file '
                                            f'{file_path}; Skipping file.')
                    continue
                except ModuleNotFoundError:
                    print('This is a test')
                if ('RG' not in bam_file.header):
                    self.vaseLogger.warning('Read group tag not present '
                                            f'in {file_path}; '
                                            'Skipping file.')
                    continue
                if (len(bam_file.header['RG']) == 0):
                    self.vaseLogger.warning(f'Read group empty in '
                                            f'{file_path}; Skipping file.')
                    continue
                if ('SM' not in bam_file.header['RG'][0]):
                    self.vaseLogger.warning('Sample group tag not present '
                                            f'in {file_path}; '
                                            'Skipping file.')
                    continue
                bam_id = bam_file.header['RG'][0]['SM']
                self.bam_map[bam_id] = file_path
                bam_file.close()
            self.vaseLogger.info('Finished scanning BAM files.')
            return self.bam_map

    def get_common_ids(self):
        self.common_ids = [vcf_id for vcf_id in self.vcf_map
                           if vcf_id in self.bam_map]
        return self.common_ids

# =============================================================================
# Example usage:
# my_map_obj = file_mapper()
# my_map_obj.file_scan(['path1', 'path2', 'path3'], 'bam')
# my_map_obj.file_scan(['pathA', 'pathB', 'pathC'], 'vcf')
# my_map_obj.get_common_ids()
# =============================================================================
