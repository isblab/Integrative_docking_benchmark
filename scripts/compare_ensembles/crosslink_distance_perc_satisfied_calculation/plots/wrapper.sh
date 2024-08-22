#!/bin/bash

python plot_crosslink_distance.py --lt5
python plot_crosslink_distance.py --b6_10
python plot_crosslink_distance.py --mt10
python plot_crosslink_distance.py --edc
python plot_crosslink_distance.py --exp
python plot_crosslink_distance.py --sel

python plot_crosslink_perc_satisfied.py --lt5
python plot_crosslink_perc_satisfied.py --b6_10
python plot_crosslink_perc_satisfied.py --mt10
python plot_crosslink_perc_satisfied.py --edc
python plot_crosslink_perc_satisfied.py --exp
python plot_crosslink_perc_satisfied.py --sel

# python calc_and_plot_native_model_xlink_dist.py summary
# python calc_and_plot_native_model_xlink_dist.py complexwise
