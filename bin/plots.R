library(rio)
setwd("~/code/EntropyStableFD/bin")
data = import('burger_1_3200.txt')
library(ggplot2)
library(ggthemes)
ggplot(data=data, aes(x=x, y=MS)) +
  geom_line()+geom_point()
library(reshape2)
melted <- melt(data, id="x")



library(rbokeh)
figure(data= melted) %>%
  ly_lines(x, value, color = variable) %>%
  ly_lines(x, REF, data = ref, legend = "REF")

library(plotly)
plot_ly(melted, x = ~x, y = ~value, color = ~variable, mode='lines')

ref = import('burger_1_reference.txt')

ggplot() +
  geom_line(data=melted, aes(x=x, y=value, color=variable)) +
  geom_line(data=ref, aes(x=x, y=REF)) +
  scale_color_ptol("cyl") +
  theme_minimal()

#Errores
errores = import('burger_1_errors.txt')
errores
Order = log2(errores[1:4,]/errores[2:5,])
Order

  errores = import('burger_2_errors.txt')
errores
Order = log2(errores[1:4,]/errores[2:5,])
Order


# Otro calculo de errores
dxr = ref$x[2] - ref$x[1]
dx = data$x[2] - data$x[1]
jj = round((data$x - ref$x[1])/dxr) + 1
xs = data.frame(x1 = data$x, x2 = ref$x[jj])
sum(abs(data$MS - ref$REF[jj]))
sum(dx*abs(data$ESC - ref$REF[jj]))
a = data$ESC
b = ref$REF[jj]
c = sum(dx*abs(a[1:800] - b[1:800]))
d = sum(dx*abs(a[2400:3200] - b[2400:3200]))

sum(dx*abs(a[1400:1800] - b[1400:1800]))

