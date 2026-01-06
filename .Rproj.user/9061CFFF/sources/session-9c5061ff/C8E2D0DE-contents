data <- read.table("cancers.txt", sep="", header=TRUE, stringsAsFactors=FALSE)
rownames(data) <- data$Type
x <- data[, -1]
x <- as.data.frame(lapply(x, as.numeric))

str(x)
boxplot(x, main="Boxplots des marqueurs de cancers")
summary(t(x))

matCors <- cor(x)
matCors
symnum(matCors, abbr.colnames=FALSE)

library(corrplot)
corrplot(matCors, type="upper", order="hclust", tl.col="black", tl.srt=45)

Xsc <- scale(x, scale=FALSE)
boxplot(Xsc, main="Boxplot des marqueurs centrés")

Sigma <- t(Xsc) %*% Xsc / nrow(Xsc)
ACP <- eigen(Sigma)

ACP$values
plot(ACP$values, type="b", xlab="Composante", ylab="Valeur propre", main="Scree plot")

inertie <- cumsum(ACP$values)/sum(ACP$values)
pourcinertie <- inertie*100
pourcinertie
plot(ACP$values/sum(ACP$values), type="b", xlab="Composante", ylab="Proportion variance expliquée")

CP <- as.matrix(Xsc) %*% ACP$vectors
rownames(CP) <- rownames(x)
colnames(CP) <- paste0("PC", 1:ncol(CP))
CP

cos2 <- rowSums(CP^2)/sum(ACP$values)
cos2

n <- nrow(Xsc)
contrib <- CP^2 / (n * matrix(ACP$values, nrow=n, ncol=ncol(CP), byrow=TRUE))
contrib

plot(CP[,1:2], pch=2, cex=cos2*5, xlab="PC1", ylab="PC2", main="Projection des cancers sur PC1-PC2")
text(CP[,1:2], labels=rownames(CP), pos=3)

Rho <- diag(1/sqrt(diag(Sigma))) %*% ACP$vectors %*% diag(sqrt(ACP$values))
rownames(Rho) <- colnames(x)
colnames(Rho) <- paste0("PC", 1:ncol(Rho))
Rho

