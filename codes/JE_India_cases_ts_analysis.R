#JE incidence timeseries data from India
#Loading library:
library(xlsx)
library(cluster)
library(ggplot2)

#Loading data:
JE_India_inc = read.xlsx("D:/OUCRU/Hannah/JEV/JEV_model/data_JE_raw/cases_sero_data/India_JE_case_timeseries.xlsx",1)
#data from: http://www.census2011.co.in/states.php
sub_India_pop_2011 = data.frame(matrix(cbind(as.character(JE_India_inc$Affected.States[-nrow(JE_India_inc)]),
                                       c(84580777,1383727,31205576, 104099452, 16787941, 1458545,25351462,
                                         32988134, 61095297, 33406061, 112374333, 2855794, 2966889, 1978502,
                                         41974218, 27743338, 72147030, 35286757, 3673917, 99812341, 10086292, 91276115)),
                                       nrow = 22, ncol = 2))
ts_only_data = JE_India_inc[-nrow(JE_India_inc),3:ncol(JE_India_inc)]

m.JE_India_inc = reshape2::melt(JE_India_inc[-nrow(JE_India_inc),-1], "Affected.States")
ggplot(m.JE_India_inc, aes(x = variable, y = value, group = Affected.States))+
  geom_path(aes(color = Affected.States))+
  annotate("text", x = "X2013", y = m.JE_India_inc$value[m.JE_India_inc$variable == "X2013"], 
           label = JE_India_inc$Affected.States[-nrow(JE_India_inc)])
#ts_only_data = ts_only_data/as.numeric(as.character(sub_India_pop_2011$X2)) #standardized data from pop country
mean_case = rowMeans(ts_only_data, na.rm = T)
median_case = sapply(data.frame(t(ts_only_data)), FUN = f <- function(x) median(x, na.rm = T))
std_case = sapply(data.frame(t(ts_only_data)), FUN = f <- function(x) sd(x, na.rm = T))
min_case = sapply(data.frame(t(ts_only_data)), FUN = f <- function(x) min(x, na.rm = T))
max_case = sapply(data.frame(t(ts_only_data)), FUN = f <- function(x) max(x, na.rm = T))

clus_data = data.frame(mean_case, median_case, std_case, min_case, max_case)
colnames(clus_data) = c("mean", "median", "std", "min", "max")
rownames(clus_data) = JE_India_inc$Affected.States[-23]

#Cluster analysis:

#Partitioning method:
kmeans(clus_data,4)
pam(clus_data, 4)

#Hierarchical method:
plot(agnes(clus_data))
plot(diana(clus_data))
