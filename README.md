# getmor_segmentation

## Dependencies and compilation

## Procedure for segmentataion of rare cels

## Programs and command line options

#### ck_segment
```
ck_segment <in_dir> <out_dir> <sample_name>
           <start_offset> <num_frames> <channel_start>
           <blur_iterations> <cell_frac>
           <high_th_base_quantile> <high_th_ratio>
           <verbose>
```
* in_dir: path to input directory
* out_dir: path to output directory
* sample_name: name of sample for output file names
* start_offset: start offset of the first frame to process
* num_frames: number of frames to process
* channel_start:
* blur_iterations: number of iterations for edge preserving Gaussian
smoothing
* cell_frac
* high_th_base_quantile
* high_th_ratio
* verbose


#### feature_extraction
```
feature_extraction <frame_dir> <label_dir> <sample_name>
                   <start_offset> <num_frames> <channel_start>
```
* frame_dir
* label_dir
* sample_name
* start_offset
* num_frames
* channel_start

#### filter_features
```
filter_features <in_features_file> <out_features_file>
                <col(1-based),min,max>*
```
* in_features_file
* out_features_file
* col(1-based),min,max

#### cell_image
```
cell_image <frame_dir> <out_dir> <cell_list_file>
           <ch_use_vector> <ch_start_offset> <ch_gain>
           <out_frame_size>
```
* frame_dir
* out_dir
* cell_list_file
* ch_use_vector
* ch_start_offset
* ch_gain
* out_frame_size



#### single_channel_mask
```
single_channel_mask <in_dir> <out_dir> <sample_name>
           <start_offset> <num_frames>
           <blur_iterations> <cell_frac> <low_th_offset>
           <high_th_base_quantile> <high_th_ratio>
           <verbose>
```
* in_dir
* out_dir
* sample_name
* start_offset
* num_frames
* blur_iterations
* cell_frac
* low_th_offset
* high_th_base_quantile
* high_th_ratio
* verbose


#### merge_and_watershed
```
merge_and_watershed <in_dir> <out_dir> <start_offset> <channel_start>
                    <num_frames> <watershed_level>
                    <verbose>
```
* in_dir
* out_dir
* start_offset
* channel_start
* num_frames
* watershed_level
* verbose


#### gain
```
gain: <in_dir> <out_dir> <start_offset> <num_frames>
      <gain>
```
* in_dir
* out_dir
* start_offset
* num_frames
* gain

## Citation

## Contact information
