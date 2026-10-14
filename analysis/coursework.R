
# Soient x: la dose, n; le nombre de scarabées, et y: le nombre de scarabées touchés

#Données et  Variables globales :

x <- c(1.6907, 1.7242, 1.7552, 1.7842, 1.8113, 1.8369, 1.8610, 1.8839)
n <- c(59, 60, 62, 56, 63, 59, 62, 60)
y <- c(6, 13, 18, 28, 52, 53, 61, 60)


p_obs <- y / n # Proportions observées de scarabées tués à chaque dose

X <- cbind(1, x) # Construction de la matrice des xi: une seule fois pour tous les modèles

tol <- 1e-6    #  # Valeur de tolérance de convergence (pour l'algorithme Newton - Raphson ): pour les 3 modèles

max_iter <- 100

## Le graphique de la proportion observée: pi=Yi/ni en fonction de la dose Xi 

plot(x, p_obs, pch = 19, col = "blue",
     main = "Proportion de scarabées tués selon la dose",
     xlab = "Dose", ylab = "Proportion tuée", ylim = c(0,1))

#=====================================================
# Commentaires sur l'allure de la relation dose-reponse:
#======================================================

# La proportion de scarabées tués augmente avec la dose d'insecticide.
# Pour les plus faibles doses, la mortalité reste faible, puis elle croît rapidement entre 1.78 et 1.83
# avant d'atteindre presque 100% pour les plus fortes doses(>1.86).

#-----------------------------------------------------
#============= Fonctions communes aux 3 modeles  ====================:
#-----------------------------------------------------

deviances <- function(y, n, p_pred) {
  2 * sum(
    ifelse(y == 0, 0, y * log(y / (n * p_pred))) +
      ifelse(y == n, 0, (n - y) * log((n - y) / (n * (1 - p_pred))))
  )
}
chi2_pearson <- function(y, n, p_pred) {
  sum((y - n * p_pred)^2 / (n * p_pred * (1 - p_pred)))
}
#=============Logit=====================
#) Implémentation de  l’algorithme de Newton-Raphson jusqu’à convergence.
beta_logit=matrix(c(0,0), ncol=1)
iter <- 1
diff <- Inf   # initialiser la différence
while(diff > tol && iter < max_iter){
  eta <- X %*% beta_logit
  p   <- 1/(1+exp(-eta))
  mu  <- n * p
  U   <- t(X) %*% (y - mu)                           # score
  W   <- diag(as.numeric(n * p * (1-p)))
  I   <- t(X) %*% W %*% X                            # Fisher
  beta_new <- beta_logit + solve(I, U)                     # mise à jour
  diff <- max(abs(beta_new - beta_logit))
  beta_logit <- beta_new; iter <- iter + 1
}
b_logit <- as.numeric(beta_logit)

cat("Estimateurs du maximum de vraisemblance obtenus :\n",b_logit)
# -------------------------
# c) PRÉDICTIONS
# -------------------------

# i) Les valeurs prédites linéaires.
eta_logit <- X %*% beta_logit
eta_logit

#ii) Les proportions prédites.
p_logit <- 1 / (1 + exp(-eta_logit))
p_logit
# iii) Les effectifs prédits pour chaque dose.
Y_logit <- n * p_logit
Y_logit

# -------------------------
# d) Visualisation des résultats :
# -------------------------
plot(x, p_obs, pch = 19, col = "blue",
     main = "Proportions observées et prédites (logit)",
     xlab = "Dose ",
     ylab = "Proportion de scarabées tués",
     ylim = c(0, 1))

# Courbe prédite logistique
curve(1 / (1 + exp(-(beta_logit[1] + beta_logit[2] * x))),
      add = TRUE, col = "red", lwd = 2)

points(x, p_logit, col = "red", pch=17)

legend("topleft",
       legend = c("Observées", "Prédites (modèle logit)"),
       col = c("blue", "red"), pch = c(19, 17), lwd = c(NA, 2))

# Commentaires sur la qualité de l’ajustement:

#Le modèle logit arrive globalement à bien suivre la forme des données. La probabilité de mortalité augmente
#bien avec la dose, et la courbe ajustée passe très proche des points, surtout pour les doses moyennes et
#élevées. Aux faibles doses, on voit que le modèle (les points en bleu ) a tendance à prédire un peu trop de mortalité, et autour
#du milieu de la plage de doses il sous-estime légèrement. Mais ces écarts restent modestes. Dans l’ensemble,
#l’ajustement est satisfaisant et la tendance est capturée correctement.


# -------------------------
# e) Diagnostics :
# -------------------------
# La deviance D pour ce modèle
D_logit <- deviances(y, n, p_logit)
D_logit

# La statistique du χ2 de Pearson
X2_logit <- chi2_pearson(y, n, p_logit)
X2_logit


#------------------------------------------------------------------------------
#=============Probit=====================
#------------------------------------------------------------------------------
beta_probit=matrix(c(0,0), ncol=1)
iter <- 1
diff <- Inf   # initialiser la différence
while(diff > tol && iter < max_iter){
  eta <- X %*% beta_probit
  p   <- pnorm(eta)
  phi <- dnorm(eta)
  mu  <- n * p
  # score et Fisher
  U   <- t(X) %*% ((y - mu) * (phi / (p*(1-p))))
  W   <- diag(as.numeric(n * (phi^2) / (p*(1-p))))
  I   <- t(X) %*% W %*% X
  beta_new <- beta_probit + solve(I, U)
  diff <- max(abs(beta_new - beta_probit))
  beta_probit <- beta_new; iter <- iter + 1
}
b_probit <- as.numeric(beta_probit)
cat("Estimateurs du maximum de vraisemblance obtenus :\n",b_probit)

# -------------------------
# c) PRÉDICTIONS
# -------------------------

# (i) Valeurs linéaires prédites
eta_probit <- X %*% beta_probit
eta_probit

#(ii) Proportions prédites 
p_probit <-pnorm(eta_probit)
p_probit
# (iii) Effectifs prédits Ŷ_i = n_i * p̂_i
Y_probit <- n * p_probit 
Y_probit

# -------------------------
# d) Visualisation des résultats :
# -------------------------
plot(x, p_obs, pch = 19, col = "blue",
     xlab = "Dose ", ylab = "Proportion tuée",
     main = "Modèle probit : Observations vs Prédictions",
     ylim = c(0,1))

curve(pnorm(beta_probit[1] + beta_probit[2]*x), add = TRUE, col = "red", lwd = 2)
points(x, p_probit, col = "red", pch = 17)

legend("topleft", legend = c("Observées", "Prédites (probit)"),
       col = c("blue", "red"), pch = c(19,17), lwd = c(NA,2))


# Commentaires sur la qualité de l’ajustement:

#L'ajustement (la courbe en rouge ) suit très bien la tendance générale des observations: 
#Comme pour le modele logit, on observe une legère surestimation de la mortalité par le modèle (les points en bleu) aux faibles doses, et une sous estimation aux doses intermediaires.
#Pour les fortes doses, l'ajustement est très bon: les points et la courbe se confondent presque .


# -------------------------
# e) Diagnostics :
# -------------------------
# La deviance D pour ce modèle
D_probit <- deviances(y, n, p_probit)
D_probit

# La statistique du χ2 de Pearson
X2_probit <- chi2_pearson(y, n, p_probit)
X2_probit


#------------------------------------------------------------------------------
#=============C-log-log=====================
#------------------------------------------------------------------------------
beta_log=matrix(c(0,0), ncol=1)
iter <- 1
diff <- Inf   # initialiser la différence
while(diff > tol && iter < max_iter){
  eta <- X %*% beta_log
  p   <- 1 - exp(-exp(eta))                 # p(eta)
  dp  <- exp(eta - exp(eta))                # dp/d(eta)
  mu  <- n * p
  
  # Score et Fisher
  U <- t(X) %*% ((y - mu) * (dp / (p * (1 - p))))
  W <- diag(as.numeric(n * (dp^2) / (p * (1 - p))))
  I <- t(X) %*% W %*% X
  
  beta_new <- beta_log + solve(I, U)
  diff <- max(abs(beta_new - beta_log))
  beta_log <- beta_new; iter <- iter + 1
}
b_log <- as.numeric(beta_log)
cat("Estimateurs du maximum de vraisemblance obtenus :\n",b_log)

# -------------------------
# c) PRÉDICTIONS
# -------------------------

# (i) Valeurs linéaires prédites
eta_log <- X %*% beta_log
eta_log

#(ii) Proportions prédites p̂_i
p_log <-1 - exp(-exp(eta_log))
p_log
# (iii) Effectifs prédits Ŷ_i = n_i * p̂_i
Y_log <- n * p_log
Y_log

# -------------------------
# d) Visualisation des résultats :
# -------------------------
plot(x, p_obs, pch = 19, col = "blue",
     xlab = "Dose ", ylab = "Proportion tuée",
     main = "Modèle c-log-log : Observées vs Prédites", ylim = c(0,1))
curve(1 - exp(-exp(beta_log[1] + beta_log[2]*x)), add = TRUE, col = "red", lwd = 2)
points(x, p_log, col = "red", pch = 17)
legend("topleft", legend = c("Observées", "Prédites (c-log-log)"),
       col = c("blue", "red"), pch = c(19, 17), lwd = c(NA, 2))


# Commentaires sur la qualité de l’ajustement:
#Le modèle complémentaire log-log suit lui aussi très bien les données. Il fournit un ajustement satisfaisant
#et cohérent avec les données observées. Les écarts sont mineurs et ne remettent pas en cause la qualité du modèle


# -------------------------
# e) Diagnostics :
# -------------------------
# La deviance D pour ce modèle
D_log <- deviances(y, n, p_log)
D_log

# La statistique du χ2 de Pearson
X2_log <- chi2_pearson(y, n, p_log)
X2_log



#---------------------------------------------------------
#=========================Comparaison des modeles: ======
#-------------------------------------------------------

res <- data.frame(
  Modele = c("Logit", "Probit", "C-log-log"),
  beta0  = c(b_logit[1],b_probit[1],b_log[1]),
  beta1  = c(b_logit[2],b_probit[2],b_log[2]),
  
  Deviance = c(D_logit,D_probit,D_log),    # valeurs calculées D
  Chi2_Pearson = c(X2_logit,X2_probit,X2_log)     # valeurs calculées X2
)
res

# Comparaison des estimateurs des paramètres obtenus pour chaque modèle:
#Les trois modèles (logit, probit et c-log-log) donnent tous un coefficient beta1_chapeau > 0, ce qui confirme que la probabilité de mortalité augmente avec la dose.
#Les modèles sont donc cohérents sur l’effet principal.

#Les valeurs diffèrent cependant selon le modèle :(probit,c-log-log ou logit). 

#Les valeurs des coefficients changent selon la fonction de lien, car elles n’utilisent pas la même échelle. Malgré cela, les trois modèles indiquent tous que la dose augmente la probabilité de mortalité.

#Sélection du modèle :
# Proposition de critère de sélection basé sur la deviance:

#Comme nos trois modèles ont le même nombre de paramètres (beta0, beta1) et sont ajustés sur les mêmes données, on peut comparer directement leurs déviances résiduelles (plus petit = meilleur ajustement)

# Avantages et inconvénients de chaque modèle:

#### Modèle logistique :

#- Avantages :
# Interprétation intuitive : les coefficients peuvent être transformés en rapport de cotes, ce qui facilite la compréhension des effets.
#- Relativement robuste aux valeurs extrêmes dans les données, comparé à un modèle linéaire classique.

# - Inconvénients :
#- Nécessite des échantillons d’au moins 10 à 20 observations pour obtenir des estimations stables et fiables.
#- L’hypothèse que les erreurs suivent une loi logistique peut être restrictive si la vraie distribution est très asymétrique.
#- Le modèle suppose une relation symétrique entre les variables explicatives et la probabilité de l’événement.


# Modèle probit :


# - Avantages :
#- Relativement robuste aux valeurs aberrantes dans les données.
#- Adapté à la modélisation de variables dépendantes discrètes.

#- Inconvénients :
#- Nécessite des échantillons de grande taille pour obtenir des estimations stables et fiables.
#- Les coefficients ne s’interprètent pas directement comme des variations marginales de la probabilité, ce qui peut rendre leur interprétation complexe.
#- L’hypothèse de distribution normale des erreurs peut ne pas être vérifiée dans certaines situations réelles, ce qui peut biaiser les résultats.


#### Modèle c-log-log :

#### - Avantages :
#- Asymétrique donc mieux adapté aux situations où la probabilité d’un événement est très faible ou très élevée.

#### - Inconvénients :
#- Suppose une distribution asymétrique des erreurs, ce qui peut ne pas convenir à toutes les données.
#- Si la relation est symétrique entre probabilité et prédicteurs, le modèle logistique ou le modèle probit peuvent mieux convenir.


# Conclusion :
##  Le  modèle qui semble le plus adapté aux données ?
#-Sur la base de la déviance résiduelle, **le modèle complémentaire log-log** obtient la valeur la plus faible parmi les trois, ce qui indique un ajustement numérique légèrement meilleur.

##  - Justification du  choix du modèle retenu en s’appuyant sur les résultats numériques et graphiques:

#- Les trois modèles décrivent correctement la relation dose–réponse et fournissent des estimateurs similaires pour l’effet de la dose (beta1_chapeau > 0). Visuellement, les ajustements sont proches et tous reproduisent bien la montée de la mortalité avec la dose.

#- D’un point de vue numérique, le modèle complémentaire log-log présente la déviance résiduelle la plus faible et la plus petite statistique de Pearson. Cela suggère un ajustement légèrement meilleur. 
#- Graphiquement, la courbe prédite par le modele c-log-log suit bien la forme observée : mortalité faible à faible dose, puis forte augmentation à forte dose;

#-En conclusion, si l’on doit sélectionner un modèle selon la déviance, **le modèle complémentaire log-log peut être retenu **. Mais les trois modèles sont cohérents et mènent à la même interprétation.