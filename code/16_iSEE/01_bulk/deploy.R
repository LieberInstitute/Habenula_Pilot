library("rsconnect")
library("here")
options(repos = BiocManager::repositories())
rsconnect::deployApp(
    appDir = here("code", "16_iSEE", "01_bulk"),
    appFiles = c("app.R", "rse_gene_Habenula_Pilot.rda"),
    appName = "habenulaPilot_bulk",
    account = "libd",
    server = "shinyapps.io"
)
