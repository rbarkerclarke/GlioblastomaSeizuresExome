caris_main = readxl::read_xlsx("~/CarisData/CARIS Data_10-27-23 with Clinical Info_deidentified.xlsx", sheet = 'Clinical Info')
caris_main$seizure_any = caris_main$`Seizure presentation (0=No, 1=Yes)` + caris_main$`Seizure presentation ONLY after surgery`
caris_main = caris_main |> dplyr::filter(`Deidentified code` %in% all$Deidentified.code)
all = readxl::read_xlsx("~/CarisData/Raw CARIS Data_8-18-24_deidentified.xlsx", sheet = 'Mutations_P_LP_VUS')
egfrs = read.csv(file = "~/CarisData/EGFR_CARIS.txt", header = FALSE)
strsplit()
library(stringr)
library(tidyr)
all <- all |> dplyr::filter(value!="NA")
  
PTEN <- all |> dplyr::filter(test=="PTEN") |> dplyr::filter(grepl("\\|", value))
#PTEN <- PTEN[grepl("|", PTEN$value), ]
PTEN <- PTEN |> separate_wider_delim(value, delim = "|", names = c("foo", "bar"))

tbl = table(PTEN$foo)
result <- paste(names(tbl), tbl, sep = "@")
paste(result, collapse=" ")


NF1 <- all |> dplyr::filter(test=="NF1") |> dplyr::filter(grepl("\\|", value))
#PTEN <- PTEN[grepl("|", PTEN$value), ]
NF1 <- NF1 |> separate_wider_delim(value, delim = "|", names = c("foo", "bar"))

tbl = table(NF1$foo)
result <- paste(names(tbl), tbl, sep = "@")
paste(result, collapse=" ")

EGFR <- all |> dplyr::filter(test=="EGFR") |> dplyr::filter(grepl("\\|", value))
#PTEN <- PTEN[grepl("|", PTEN$value), ]

TP53 <- all |> dplyr::filter(test=="TP53") |> dplyr::filter(grepl("\\|", value))
#PTEN <- PTEN[grepl("|", PTEN$value), ]
TP53 <- TP53 |> separate_wider_delim(value, delim = "|", names = c("foo", "bar"))

tbl = table(TP53$foo)
result <- paste(names(tbl), tbl, sep = "@")
paste(result, collapse=" ")


PTEN <- all |> dplyr::filter(test=="PTEN") |> dplyr::filter(grepl("\\|", value))
#PTEN <- PTEN[grepl("|", PTEN$value), ]
PTEN <- PTEN |> separate_wider_delim(value, delim = "|", names = c("foo", "bar"))

tbl = table(PTEN$foo)
result <- paste(names(tbl), tbl, sep = "@")
paste(result, collapse=" ")

SETD2 <- all |> dplyr::filter(test=="SETD2") |> dplyr::filter(grepl("\\|", value))
#PTEN <- PTEN[grepl("|", PTEN$value), ]
SETD2 <- SETD2 |> separate_wider_delim(value, delim = "|", names = c("foo", "bar"))

tbl = table(SETD2$foo)
result <- paste(names(tbl), tbl, sep = "@")
paste(result, collapse=" ")

PIK3CA <- all |> dplyr::filter(test=="PIK3CA") |> dplyr::filter(grepl("\\|", value))
#PTEN <- PTEN[grepl("|", PTEN$value), ]
PIK3CA <- PIK3CA |> separate_wider_delim(value, delim = "|", names = c("foo", "bar"))

tbl = table(PIK3CA$foo)
result <- paste(names(tbl), tbl, sep = "@")
paste(result, collapse=" ")


tbl = table(egfrs$foo)
result <- paste(names(tbl), tbl, sep = "@")

# Print result
print(result)
paste(result, collapse=" ")


sort(table(all$test))


## THOUGHTS -> BINOMIAL 
## SEIZURE Y/N SIGNIFICANT
PIK3CA_res = caris_main |> dplyr::filter(`Deidentified code` %in% unique(PIK3CA$Deidentified.code)) |> dplyr::select(seizure_any)
p_success = sum(caris_main$seizure_any, na.rm = TRUE)/length(caris_main$seizure_any)
  
binom.test(x = sum(PIK3CA_res$seizure_any), n=length(PIK3CA_res$seizure_any), p=p_success)
p_success
