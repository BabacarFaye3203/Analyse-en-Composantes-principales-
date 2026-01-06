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



#exo 2

# 1️⃣ Lecture et préparation des données
protein <- read.csv("protein.csv", header=TRUE, stringsAsFactors=FALSE)
rownames(protein) <- protein$Country     # mettre le nom des pays en ligne
x <- protein[, -1]                       # retirer la colonne Country
x <- as.data.frame(lapply(x, as.numeric)) # convertir en numérique

# Vérification
str(x)
boxplot(x, main="Boxplots de consommation alimentaire par pays")
summary(t(x))  # résumé par pays (individus)

# 2️⃣ Matrice de corrélation
matCors <- cor(x)
matCors
library(corrplot)
corrplot(matCors, type="upper", order="hclust", tl.col="black", tl.srt=45)

# 3️⃣ ACP avec prcomp
res.pca <- prcomp(x, center=TRUE, scale.=TRUE) # centrer et standardiser

# Valeurs propres et proportion de variance expliquée
res.pca$sdev^2                    # valeurs propres
prop_var <- (res.pca$sdev^2)/sum(res.pca$sdev^2)
cumsum(prop_var)                  # inertie cumulée
plot(prop_var, type="b", main="Variance expliquée par chaque PC", xlab="PC", ylab="Proportion")

# 4️⃣ Composantes principales
CP <- res.pca$x                     # coordonnées des pays (individus)
cos2 <- rowSums(CP^2)/sum(res.pca$sdev^2)
CP[1:5,]                            # afficher les 5 premiers pays

# 5️⃣ Contributions
n <- nrow(x)
contrib <- CP^2 / (n * matrix(res.pca$sdev^2, nrow=n, ncol=ncol(CP), byrow=TRUE))
contrib[1:5,]                        # contributions pour les 5 premiers

# 6️⃣ Représentation graphique
plot(CP[,1:2], pch=16, cex=cos2*5, xlab="PC1", ylab="PC2", main="Projection des pays sur PC1-PC2")
text(CP[,1:2], labels=rownames(CP), pos=3)

# 7️⃣ Corrélation des variables avec les PC
Rho <- cor(x, CP)
Rho

# 8️⃣ ACP avec FactoMineR pour plus d’options
library(FactoMineR)
library(factoextra)

res.acp <- PCA(x, graph=FALSE)

# Valeurs propres et scree plot
eig.val <- get_eigenvalue(res.acp)
eig.val
fviz_eig(res.acp, addlabels = TRUE)

# Variables
var <- get_pca_var(res.acp)
head(var$coord)      # coordonnées
head(var$cos2)       # qualité de représentation
head(var$contrib)    # contributions

# Graphiques
fviz_pca_var(res.acp, col.var="cos2", gradient.cols=c("#00AFBB","#E7B800","#FC4E07"), repel=TRUE)
fviz_contrib(res.acp, choice="var", axes=1:2, top=10)
fviz_pca_biplot(res.acp, repel=TRUE, col.var="#2E9FDF", col.ind="#696969")

