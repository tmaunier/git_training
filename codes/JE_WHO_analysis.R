#JE incidence data from WHO classification:
#Loading library:
library(xlsx)
library(cluster)

#Loading data:
JE_WHO_inc = read.xlsx("D:/OUCRU/Hannah/JEV/JEV_model/data_JE_clean/JE_inci_sero/JE_incidence_time_series_WHO.xlsx",1)
sel_countries = c("BGD", "BTN", "KHM", "PRK", "IND", "IDN", "LAO", "MMR", "NPL", "PNG", "LKA", "TLS", "VNM")

Inc_sel = JE_WHO_inc[match(sel_countries, JE_WHO_inc$ISO_code),]
ts_only_data = Inc_sel[,5:ncol(Inc_sel)]

mean_case = rowMeans(ts_only_data, na.rm = T)
median_case = sapply(data.frame(t(ts_only_data)), FUN = f <- function(x) median(x, na.rm = T))
std_case = sapply(data.frame(t(ts_only_data)), FUN = f <- function(x) sd(x, na.rm = T))
min_case = sapply(data.frame(t(ts_only_data)), FUN = f <- function(x) min(x, na.rm = T))
max_case = sapply(data.frame(t(ts_only_data)), FUN = f <- function(x) max(x, na.rm = T))

clus_data = data.frame(mean_case, median_case, std_case, min_case, max_case)
colnames(clus_data) = c("mean", "median", "std", "min", "max")
rownames(clus_data) = Inc_sel$ISO_code

#Cluster analysis:

#Partitioning method:
kmeans(clus_data, 4)
pam(clus_data, 4)

#Hierarchical method:
plot(agnes(clus_data))
plot(diana(clus_data))
