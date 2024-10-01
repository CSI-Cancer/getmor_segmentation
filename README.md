# getmor_segmentation

GeTMoR (GEnomic, Transcriptomic, and MOrphological profiling of Rare cells)
is a protocol for or detecting rare cancer related cells to
simultaneously image and profile the genome and transcriptome from single
rare cells. This repo contains the code that we used for image analysis
of scanned Immunofluorescence slides.

## Dependencies and compilation
1. We make extensive use of the ITK toolkit for image analysis.
Download and install the
[ITK toolkit](https://itk.org/):

2. Set the `ITK_DIR` variable in `src/CMakeLists.txt`

3. Create a build director:
```
mkdir build; cd build
```

4. Create makefiles using:
```
cmake <path_to_rerp>/getmor_segmentation/src/
```

5. Compile the code with:
```
make
```


## Procedure for GeTMoR rare cell detection
The parameters for the programs below are based on the assumption that
the cells in the cytokeratin (CK) channel are rare and so at most 3 cells
are present in a frame. If this assumption is not met, the parameters
need to be changed accordingly to achieve high sensitivity of detection.

It is assumed that the scanned image file names are of the format
`Tile\%06d.tif`, that 2304 images are generated per channel, and that
the different channel images are named sequentially after each other. If
differing number of frames are imaged, adjust the frame offsets in the
commands accordingly.

* Run the command below in a cmd prompt window after 20\% of the
  frames are scanned on the CK channel to segment the cell:
```
 ck_segment <scan_output_dir> <segmentation_dir> <sample_name>
 1 2304 2305 5 0.995 0.3 3 1
 ```
 where `scan_output_dir` is the path to the directory
 containing scanned images, `segmentation_dir` is the output
 directory to store the segmented images, and `sample_name` is an
 identifier for the sample.  


 * Run the command to extract cell features:
 ```
 $ feature_extraction <scan_output_dir> <segmentation_dir>
 <sample_name> 1 2304 <channel_start>
 ```
 where `channel_start` is a comma separated list of channel start
 offsets to user for feature extraction. This should be set to
 ``"1, 2305"`` when scanned on DAPI and CK channel and to
 ``"2305"`` when scanned only on the CK channel.

 * Run the command to filter cells based on feature values:
 ```
 $ filter_features <sample_name>_feature_vec.txt
 <sample_name>_feature_vec_filt.txt <col(1-based),min,max>*
 ```
 where ``<col(1-based),min,max>*`` specifies the column in the
 feature vector file and the minimum and maximum values in that column to
 retain. To filter events corresponding to typical cell sizes use
 ``"5,100,10000"``. In addition, filters can also be used for DAPI
 and CK mean intensities.

 * Run the command to visualize and confirm all the segmented
 events:
 ```
 $ cell_image <scan_output_dir> <segmentation_dir>
 <sample_name>_feature_vec_filt.txt "1,1,0,0" "2305,1,0,0" 150
 ```


## Citation

## Contact information
Rishvanth K. Prabakar kaliapp@cshl.edu
