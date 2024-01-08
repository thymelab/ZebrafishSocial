#!/bin/bash
#SBATCH -p express # Partition to submit to
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 60 # Runtime in minutes
#SBATCH --mem=64000 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH -o linearmodel_hcn1_hetandwt_vs_hcn1_hom_%A_%a.out # Standard out goes to this file
#SBATCH -e hostname_%A_%a.err # Standard err goes to this filehostname
module load Anaconda3
cd output_mergingfinal
python /data/project/thymelab/behaviordata/socialbehavior/socialscripts/processmotiondata.py -t "../testlog.timestamp1.Fri, Apr 8, 2022" -e "../fulltestrun_final_01_27_2020" -c "../testlog.centroid1.Fri, Apr 8, 2022" -d "../testlog.motion1.Fri, Apr 8, 2022" -m "../hsmovieFri, Apr 8, 2022_" -g "../hcn1_hetandwt_vs_hcn1_hom_scripted_inputgenotypeids" -s "../sectionsfile" -j "../../../scripts/behavior_published2020/prepandsupportfiles/PlotParameters" -n 20 -r "../rois_string" -graphmulti "/data/project/thymelab/behaviordata/socialbehavior/hcn1_social_plate1_Box10_4_8_2022/outputfulldata_hcn1_hetandwt_vs_hcn1_hom/,/data/project/thymelab/behaviordata/socialbehavior/hcn1_social_plate2_Box11_4_8_2022/outputfulldata_hcn1_hetandwt_vs_hcn1_hom/,/data/project/thymelab/behaviordata/socialbehavior/hcn1_social_plate3_Box12_4_8_2022/outputfulldata_hcn1_hetandwt_vs_hcn1_hom/"
