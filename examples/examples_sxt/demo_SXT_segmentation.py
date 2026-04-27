#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Examples for calling SXT cell/mitochondria/insulin vesicle segmentation script
"""
#%%
# cell_segmentation

import os
import sys
import argparse
from datetime import datetime
import logging
from multiprocessing import Pool
import warnings

# Filter MRC file warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, module="mrcfile")
warnings.filterwarnings("ignore", message=".*Unrecognised machine stamp.*")
warnings.filterwarnings("ignore", message=".*Map ID string not found.*")
warnings.filterwarnings("ignore", message=".*data block cannot be read.*")

from ipa.processing.segmentation import run_cell_segmentation
from ipa.processing.segmentation import run_mito_segmentation
from ipa.processing.segmentation import run_sphere_like_organelle_segmentation
from ipa.processing.segmentation import train_net, UNet
import torch

from ipa.data_loader import UniversalDataLoader
from ipa.data_loader import QuickLogger
from parsers import get_args
import glob
import time
mainpath = get_args().main_path



def cell_segmentation():


    mainpath = get_args().main_path

    # Initialize logger
    log_dir = f'{mainpath}/logs'
    logger = QuickLogger("sxt_cell_segmentation", log_dir=log_dir)
    logger.step("Starting SXT Cell Segmentation Demo")

    dataid_list = ['784_5']
    image_data = {}
    logger.step(f"Processing datasets: {dataid_list}")
    

    for dataid in dataid_list:
        image_path = f'{mainpath}/data/sxt_images/Stevens_pancreatic_INS_1E_784_5_pre_rec.mrc'
        logger.file_in(image_path)
        logger.step(f"Loading data for {dataid} from: {image_path}")

        data = UniversalDataLoader.load_data(image_path)
        image_data[dataid] = data  # Extract image data part

        logger.step(f"Successfully loaded {dataid}: {image_data[dataid].shape}, dtype: {image_data[dataid].dtype}")

    # model_path = 'models'  # model form your path
    save_dir = f'{mainpath}/results/mem_nu_seg_with_val'
    logger.step(f"Save directory: {save_dir}")

    logger.step("Running cell segmentation...")
    run_cell_segmentation(
        # model_path=model_path,
        save_dir=save_dir,
        pool_processes=6,
        dataid=dataid_list,
        image_data=image_data  
    )
    
    logger.step("Cell segmentation completed successfully")


#%%
# mitochondria_segmentation


def mito_segmentation():


    mainpath = get_args().main_path

    # Initialize logger
    log_dir = f'{mainpath}/logs'
    logger = QuickLogger("sxt_mito_segmentation", log_dir=log_dir)
    logger.step("Starting SXT Mitochondria Segmentation Demo")

    dataid_list = ['784_5']
    image_data = {}
    logger.step(f"Processing datasets: {dataid_list}")

    for dataid in dataid_list:
        image_path = f'{mainpath}/data/sxt_images/Stevens_pancreatic_INS_1E_784_5_pre_rec.mrc'
        logger.file_in(image_path)
        logger.step(f"Loading data for {dataid} from: {image_path}")

        data= UniversalDataLoader.load_data(image_path)
        image_data[dataid] = data  

        logger.step(f"Successfully loaded {dataid}: {image_data[dataid].shape}, dtype: {image_data[dataid].dtype}")

    # model_path = 'models'  # model form your path
    save_dir = f'{mainpath}/results/mem_nu_seg_with_val'
    logger.step(f"Save directory: {save_dir}")

    logger.step("Running mitochondria segmentation...")
    run_mito_segmentation(
        # model_path=model_path,
        save_dir=save_dir,
        pool_processes=6,
        dataid=dataid_list,
        image_data=image_data  # Pass pre-loaded image data
    )
    
    logger.step("Mitochondria segmentation completed successfully")


def isg_segmentation():


    mainpath = get_args().main_path

    # Initialize logger
    log_dir = f'{mainpath}/logs'
    logger = QuickLogger("sxt_isg_segmentation", log_dir=log_dir)
    logger.step("Starting SXT ISG Segmentation Demo")

    dataid_list = ['784_5']
    image_data = {}
    logger.step(f"Processing datasets: {dataid_list}")

    for dataid in dataid_list:
        image_path = f'{mainpath}/data/sxt_images/Stevens_pancreatic_INS_1E_784_5_pre_rec.mrc'
        logger.file_in(image_path)
        logger.step(f"Loading data for {dataid} from: {image_path}")

        data= UniversalDataLoader.load_data(image_path)
        image_data[dataid] = data  

        logger.step(f"Successfully loaded {dataid}: {image_data[dataid].shape}, dtype: {image_data[dataid].dtype}")

    # model_path = 'models'  # model form your path
    output_dir = f'{mainpath}/results/isg_semantic_mask'
    logger.step(f"Output directory: {output_dir}")

    logger.step("Running ISG segmentation...")
    run_sphere_like_organelle_segmentation(
        save_dir=output_dir,
        pool_processes=1,
        dataid=dataid_list,
        image_data=image_data,
        # model_path=model_path,
    )
    
    logger.step("ISG segmentation completed successfully")



if __name__ == "__main__":
    cell_segmentation()
    mito_segmentation()
    isg_segmentation()


