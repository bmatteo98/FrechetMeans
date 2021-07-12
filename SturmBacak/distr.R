
df <- data.frame(
  sex=factor(rep(c("F", "M"), each=200)),
  weight=round(c(rgamma(200, 5, 0.5), rgamma(200, 2, 0.3)))
)
head(df)

library(ggplot2)
# Basic histogram
ggplot(df, aes(x=weight))+
  geom_histogram(color="darkblue", fill="lightblue")+
  geom_vline(aes(xintercept=mean(weight)), color="red",
             linetype="dashed", size = 1.3)+
  geom_vline(aes(xintercept=median(weight)), color="darkgreen",
             linetype="dashed", size = 1.3)+
  geom_vline(aes(xintercept=5.8), color="darkviolet",
           linetype="dashed", size = 1.3)
