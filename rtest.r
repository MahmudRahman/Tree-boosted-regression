Quora = read.table('quoraforR.txt',head = F)
Quora[Quora[,22] <0, 22] = 0
library(randomForest)
library(quadprog)
#library(kernlab)
#library('LowRankQP')

Quora.raw = Quora
Quora.std = Quora.raw
library(pls)
Quora.std[,c(1,2,3,4,5,11,20,21)] = stdize(Quora.raw[,c(1,2,3,4,5,11,20,21)])
Quora.std$V6 = Quora.raw$V6/max(Quora.raw$V6)

#### 10 folds CV ######
n = dim(Quora.std)[1]
p = dim(Quora.std)[2] - 1
Quora.std[,22] = as.factor(Quora.std[,22])
V = floor(n/10)
m = 250 #tree number
r1 = c(0, 0.15,  0.2, 0.25, 0.3, 0.4) #fraction of err estimating data
Amat = t(rbind(rep(1,m), diag(rep(1,m),m)))
b0 = c(1,rep(0,m))

CVerr3 = CVerr2 = CVerr = CVerr.true = CVerr.false = matrix(0,1+length(r1), 10)

for(i in 1:10){
    #10 fold CV
    test.index = c((V*(i-1)+1):min((V*i), n))
    data.test = Quora.std[test.index,-22]
    data.train = Quora.std[-test.index,-22]
    Y.test = Quora.std[test.index, 22]
    Y.train = Quora.std[-test.index, 22]

    # Random Forest
    Res = randomForest(x = data.train, y = Y.train, replace = FALSE, keep.forest = TRUE, ntree = m)

    pred.test.whole = matrix(as.numeric(predict(Res, newdata = data.test, predict.all = TRUE)$individual), ncol = m)%*%matrix(rep(1,m)/m, ncol = 1)
    err.test.whole = mean((pred.test.whole>0.5) == (Y.test==1))
    CVerr[1,i] = err.test.whole
    
    #logistic forest
    E1 = TBR(x = as.matrix(data.train), ntree = m, mtry = 11, y = Y.train, nodesize = 3, ridge = FALSE)
    E2 = TBR(x = as.matrix(data.train), ntree = m, mtry = 11, y = Y.train, nodesize = 3, ridge = TRUE)
    P1 = predict.TBR(ensamble = E1, newdata = as.matrix(data.test), level = levels(Y.test))
    P2 = predict.TBR(ensamble = E2, newdata = as.matrix(data.test), level = levels(Y.test))
    CVerr[2,i] = mean(P1$class == Y.test)
    CVerr[3,i] = mean(P2$class == Y.test)
    
    #singla logistic tree
    A1 = Tree_build(x = as.matrix(data.train), mtry = 21, nodesize = 3, y = Y.train, ridge = FALSE)
    A2 = Tree_build(x = as.matrix(data.train), mtry = 21, nodesize = 3, y = Y.train, ridge = TRUE)
    B1 = predict_tree(A1, newdata = as.matrix(data.test),level = levels(Y.test))
    B2 = predict_tree(A2, newdata = as.matrix(data.test),level = levels(Y.test))
    CVerr[4,i] = mean(B1==Y.test)
    CVerr[5,i] = mean(B2==Y.test)

    #original logistic regression
    M1 = glm(yy~., data = data.frame(yy=Y.train, data.train), family = binomial(link="logit"))
    CVerr[6,i] = mean(levels(Y.test)[(predict(M1, newdata = data.frame(data.test), type = 'response')>0.5)+1] == Y.test)
    
    print(i)
}
