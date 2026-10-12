# Train ChromBPNet models

`run_chrombpnet_models.py` contains the Parsl model-training task and its submission loop. It trains one model per bulk sample and chromosome fold using the selected pre-trained bias model.

The selected-bias mapping is retained in the task code as used in the manuscript workflow, including reuse of the B-cell D1 OSKM+A bias model for the B-cell conditions where that mapping was used.
