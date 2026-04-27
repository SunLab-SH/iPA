Instance segmentataion for insulin vesicle.   

## Workflows:  

### Step1:  
Generate crop image dataset.  
./crop/crop.py  
  
### Step2:  
Generate instance label.  
./instance_label_generation/gen_ins_label.py  
Validation results in ./instance_label_generation/validation_results    

### Step3:  
Denoise image.  
./denoise/predictn2v.py  

#### dependencies:
n2v = 0.11.1
tensorflow = 1.15

### Step4:  
Split image into small region with several insulin vesicles.  
./split_regions/split_img.py  
Validation results in ./split_regions/split_verify_poly  
./split_regions/split_verify_rle  

### Step5: 
Model predictions.  
Codes are in ./upsnet  

### Step6:
Merge prediction results back to original size image, opposite to splt in step 4.    
Threashold IoU = 0.5, if the IoU of current image to the near image larger than the threasold, the predicted region will be reserved.  
Codes are in ./merge/merge.py  


### Step7:
Uncropp image. opposite to crop in step 1. 
Codes are in ./resize2original_imagesize/feedback/


### Step8: 
3D fusion. Fuse 2d image to 3D image. The results are matched with SXT.
Codes are in ./3d_fusing. 3D.py, then tiff_3d.py.

## Instance Results
| Cell ID | Insulin vesicle |
|:---------------------:|:---------------------------:|
| 766_8 | 84.20|
| 784_5 | 87.56|
| 842_17 | 95.51||
| mean | **89.09** |



## Reference:
Instance segmentation frameworks are partially based on [UPSnet](https://github.com/charlotte12l/UPSNet)



## step00:
Convert raw image to recognized input data.
将原图的intensity转为0-255的图像，每一张2d单独进行转换。
Code in  convert_inputimg

