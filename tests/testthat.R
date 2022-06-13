
library(CoQR)

data_Assets <- readRDS(file = "../replication_CoQR/application/data/data_Assets.rds")

# 0. Set options
Date_max <-  "2011-12-05"

CoCAViaR_model_set <- c("CoCAViaR_SAV_diag", "CoCAViaR_SAV_fullA", "CoCAViaR_SAV_full",
                        "CoCAViaR_AS_pos", "CoCAViaR_AS_signs", "CoCAViaR_AS_mixed")


CoCAViaR_modelfits <- list()
for (jj in 1:length(CoCAViaR_model_set)){
  set.seed(1)
  CoCAViaR_model <- CoCAViaR_model_set[jj]
  Mest_obj <- CoQR(x=data_Assets %>% filter(Asset=="JPM", Date <= Date_max) %>% pull(NegReturn),
                   y=data_Assets %>% filter(Asset=="SP500", Date <= Date_max) %>% pull(NegReturn),
                   z=NULL, model=CoCAViaR_model, SRM="CoVaR", beta=0.95, alpha=0.95)
  CoCAViaR_modelfits[[paste(CoCAViaR_model)]] <- summary(Mest_obj)
}





CoCAViaR_modelfits$CoCAViaR_6p
CoCAViaR_modelfits$CoCAViaR_8pCrossA
CoCAViaR_modelfits$CoCAViaR_9p

CoCAViaR_modelfits$CoCAViaR_8pPos
CoCAViaR_modelfits$CoCAViaR_10pSigns
CoCAViaR_modelfits$CoCAViaR_10pSignsAbs

