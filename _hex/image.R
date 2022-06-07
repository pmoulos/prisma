library(ggplot2)
library(ggimage)

# A simple double density plot
df <- data.frame(
    His=c(rnorm(1000),rnorm(1000)+2),
    Cat=c(rep("First",1000),rep("Second",1000))
)
# The human
di <- data.frame(x=4,y=0.25,Cat="Third",image="human.png")

# The plot
gg <- ggplot(df,aes(x=His,colour=Cat,fill=Cat)) + 
    geom_density(alpha=0.3) +
    geom_vline(xintercept=-1.5,linetype="dashed",size=1,colour="grey40") +
    theme_void() +
    scale_color_manual(values=c(
        `First`="red3",
        `Second`="green3"
    )) +
    scale_fill_manual(values=c(
        `First`="red3",
        `Second`="green3"
    )) +
    geom_image(data=di,mapping=aes(x=x,y=y,image=image),colour="skyblue",
        size=0.3) +
    theme(legend.position="none")

ggsave(filename="image.png",plot=gg)

# Build with hexmaker - https://connect.thinkr.fr/hexmake/
# and the archived hex file
