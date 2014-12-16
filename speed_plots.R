library("ggplot2")
library("scales")

for_each_test <- read.csv("~/Dropbox/JBParallel/for_each_test.txt", header=FALSE)
min_test <- read.csv("~/Dropbox/JBParallel/min_test.txt", header=FALSE)
reduce_test <- read.csv("~/Dropbox/JBParallel/reduce_test.txt", header=FALSE)
sort_test <- read.csv("~/Dropbox/JBParallel/sort_test.txt", header=FALSE)

ggplot(for_each_test, aes(V1)) +  
  scale_y_continuous(limits = c(0,6)) +
  scale_x_continuous(trans=log2_trans(),limits=c(10000,100000000)) +
  geom_line(aes(y = V2, color = "STL for_each loop")) + 
  geom_line(aes(y = V3, color = "OpenMP For Loop")) +
  geom_line(aes(y = V4, color = "JBParellel for_each Loop")) + 
  labs(title ="Comparison of for_each loop performance", x = "Number  of elements in vector", y = "Time taken (seconds, logarithmic)")

ggplot(sort_test, aes(V1)) +  
  scale_y_continuous(limits = c(0,11)) +
  scale_x_continuous(trans=log2_trans(),limits=c(10000,100000000)) +
  geom_line(aes(y = V2, color = "STL sort loop")) + 
  geom_line(aes(y = V3, color = "JBParallel Sort")) +
  labs(title = "Comparison of sort performance", x = "Number  of elements in vector", y = "Time taken (seconds, logarithmic)")

ggplot(reduce_test, aes(V1)) +  
  scale_y_continuous(limits = c(0,0.2)) +
  scale_x_continuous(trans=log2_trans(),limits=c(10000,100000000)) +
  geom_line(aes(y = V2, color = "STL Reduce")) + 
  geom_line(aes(y = V3, color = "JBParallel Reduce")) +
  labs(title ="Comparison of reduce loop performance", x = "Number  of elements in vector", y = "Time taken (seconds, logarithmic)")

ggplot(min_test, aes(V1)) +  
  scale_y_continuous(limits = c(0,0.15)) +
  scale_x_continuous(trans=log2_trans(),limits=c(10000,100000000)) +
  geom_line(aes(y = V2, color = "STL Min")) + 
  geom_line(aes(y = V3, color = "JBParallel Min")) +
  labs(title = "Comparison of Min performance", x = "Number  of elements in vector", y = "Time taken (seconds, logarithmic)")


grid0normal <- read.csv("~/Dropbox/JBParallel/grid0normal.txt", header=FALSE)
ggplot(data=grid0normal, aes(x=V1, y=V3, fill=V2)) + geom_bar(stat="identity", position=position_dodge()) +
	  labs(title = "Comparison of Symplectic Euler Step in Mass_spring", x = "Grid Number", y = "Time (Seconds)", fill = "for_each implementation")
