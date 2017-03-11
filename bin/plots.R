library(rio)
data = import('output.txt')
library(ggplot2)
library(ggthemes)
ggplot(data=data, aes(x=x, y=MS)) +
  geom_line()+geom_point()
library(reshape2)
melted <- melt(data, id="x")
ggplot(data=melted, aes(x=x, y=value, color=variable)) +
  geom_line() + scale_color_ptol("cyl") +
  theme_minimal()


library(rbokeh)
figure(data= melted) %>%
  ly_lines(x, value, color = variable)

library(plotly)
plot_ly(melted, x = ~x, y = ~value, color = ~variable, mode='lines')
  