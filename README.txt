Currently: 2 RF models trained on inhibitor and activator bioactivity data to predict drug candidate binding for peroxisome proliferator-activated receptor gamm (PPARG), with a simple CLI to input a potential drug candidate's SMILES code. 

=====================================================================================

I did a more thorough write up of this project on my substack! 

Pt 1: background on pparg, data preprocessing, fixing class imbalance
https://substack.com/home/post/p-194485123

Pt 2: training, evalution, and validation (and realizing my architecture was naive lol)
https://laurenlog.substack.com/p/pparg-drug-candidate-prediction-pipeline-17b

My plan for the third iteration: 
Still 2 RF models, but this time one for active vs inactive compounds, and a second for agonists vs antagonist compounds. I also really want to build out an interactive UI for this project with molecular structure and data visualizations, post-predictive screening with ligand efficiency metrics, and a docking simulation, but that's gonna take me a miniute. 
