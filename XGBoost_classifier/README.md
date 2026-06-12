# Panel_2_XTi



Supervised cell type classifier for Panel 2 IMC data. 

## Workflow

   steinbock-derived single-cell measurements, 
   compensates marker spillover, 
   selects representative images for manual annotation, 
   consolidates gated labels into cell-level labels, 
   trains an XGBoost classifier.

## Concept

The classifier is trained on manually labelled cells from a curated set of 68 (49 main + 19 extra) images. 
These images are selected to cover the relevant cohorts, panels, patients, and rare cell populations. 
A separate set of 30 additional images is labelled mainly to test rare or difficult cell types outside the main training split.

The main data object throughout the workflow is a `SpatialCellExperiment`, 
where columns are cells and rows are IMC markers. 
Cell metadata stores image, patient, cohort, panel, object number, and final cell labels.

## Workflow

1. `01_read_data.Rmd`  
   Reads steinbock output for the three Panel 2 batches and prepares per-panel SPE objects. 
   It also loads image and mask objects used later for manual inspection and gating.

2. `02_spillover_correction.Rmd`  
   Applies spillover compensation to marker intensities. 
   The compensated objects are the basis for image selection, visualization, labeling, and classifier training.

3. `03_select_samples_classifier.Rmd`  
   Selects the initial 49 images across cohorts and panels. 
   These images form the main manually labelled set. 
   Prepares matching image/mask subsets and a UMAP/QC object to inspect whether the selected images cover the dataset reasonably.

4. `04_extra_20_30_imgs_SCE.Rmd`  
   Adds 19 more images to improve coverage of underrepresented cohorts and cell types, giving the final 68-image training pool. 
   It also selects 30 additional images that are kept separate and used later for rare-cell labelling and independent testing.

5. `05_consolidate_labels.Rmd`  
   Imports manual gate labels, maps them back to cells by `sample_id` and `ObjectNumber`, 
   Resolves duplicated or conflicting labels, removes problematic cells, 
   Writes labelled SPE objects for the 68-image and 30-image sets.

6. `06_train_classifier_68_main_30_to_test_imgs.Rmd`  
   ## Classifier Training
   Only cells with a confident manual label are used for training. 
   Cells labelled as `unlabelled` or with missing labels are excluded from the training matrix. 
   Before training, a few biologically inconsistent labels are cleaned
   E.g., low-IgG `Igg` cells, low-CD7 `NK` cells, and mregDC-like cells with strong epithelial/tumor marker signal.
   
   The feature matrix is built from marker expression values in the `exprs` assay. 
   DNA and histone channels are excluded so the classifier learns from biologically meaningful marker expression. 
   Rows in the model matrix are cells; columns are retained protein markers.
   
   Cell type labels are converted into numeric class IDs for XGBoost, 
   while the original label mapping is saved separately so predictions can be translated back into biological cell type names.
   
   Model selection is done with grouped cross-validation (cv) so cells from the same image are kept together during splitting. 
   This avoids overly optimistic performance estimates caused by training and testing on cells from the same image. 
   Several XGBoost hyperparameter combinations are compared using multiclass log loss, and the best-performing parameter set is selected.
   
   The final model is then retrained on all available labelled training cells from the 68-image set using the selected parameters and number of boosting rounds. 
   This final model is the classifier intended for downstream prediction.
   
   The labelled 30-image object, is kept separate from the main training set. 
   It is used as an additional held-out evaluation set
   Especially for rare or difficult classes such as `BnT`, `Igg`, `MDSC`, `mregDC`, `Neutrophil`, and `vCAF`.
   
   Training object: `spe_label_68.rds`
   Held-out/rare test object: `sce_label_30.rds`
   
   Classifier outputs:

  `classifier_xgboost_final.rds`: R-native serialized model copy used for convenience in some plotting/QC chunks.
  
  `classifier_xgboost_final.json`: portable XGBoost model file
  `classifier_xgboost_final_meta.rds`: metadata needed to use the model correctly, including label mapping, feature names, selected parameters, and number of rounds.
  
  The metadata file is important: predictions should always use the same marker order stored in `feature_names`, 
  and predicted numeric classes should always be decoded using the saved `label_mapping`.
  
  The JSON model + metadata are the safer portable option if the RDS model causes loading issues.
