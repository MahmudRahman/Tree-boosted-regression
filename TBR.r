#This is the primary function for building up a bunch of logistic trees
TBR <- function(x, y, ntree = 100, mtry = NULL, nodesize = 1, method = 'logit', ridge = TRUE){
    n = dim(x)[1]
    p = dim(x)[2]
    if(is.null(mtry)){
        mtry = ceiling(sqrt(p))
    }

    ensamble = list()
    for(k in 1:ntree){
        #build the tree node
        ensamble[[k]] = Tree_build(x = as.matrix(x), y = y, mtry = mtry, nodesize = nodesize, method = method, ridge = ridge)
        #print(k)
    }
    return(ensamble)
}

#This is the predict function for TBR object
predict.TBR <- function(ensamble, newdata, level){
    if(is.null(dim(newdata))){
        n = 1
    }else{
        n = dim(newdata)[1]
    }

    m = length(ensamble)
    vote = matrix(0, n, m)
    for(k in 1:m){
        vote[,k] = predict_tree(ensamble[[k]], newdata = newdata, level = level)
    }
    predicted = (vote == level[2])%*%matrix(1/m, m, 1)
    class = level[1+(predicted>0.5)]
    return(list(prob = predicted, class = class))
}


#This code builds up a single logistic tree
Tree_build <- function(x, y, mtry = NULL, nodesize = 1, method = 'logit', ridge = TRUE){
    require(glmnet)
    if(is.null(dim(x))){
        n = 1
        p = length(x)
    }else{
        n = dim(x)[1]
        p = dim(x)[2]
    }
    if(is.null(mtry)){
        mtry = ceiling(sqrt(p))
    }

    tree = list(value = NULL, feature = NULL, left = NULL, right = NULL)

    if(n <= nodesize){
        tree$value = levels(y)[1+(mean(y == levels(y)[2])>0.5)]
        return(tree)
    }else if(sum(y==levels(y)[1]) == n||sum(y==levels(y)[2]) == n){
        tree$value = levels(y)[1+(mean(y == levels(y)[2])>0.5)]
        return(tree)
    }else{
        #random select features
        feature = sample(1:p, mtry)
        while(all(diag(cov(x[,feature])) == 0)){
            feature = sample(1:p, mtry)
        }
        tree$feature = feature
        #print(feature)

        #run logistic
        #print(c(n, p))
        if(ridge == TRUE||n<=2*mtry){
            M1 = glmnet(x = x[,feature], y = y, family = "binomial", alpha = 0.5, lambda = 0.01)
            yhat = predict(M1, newx = x[,feature], type = 'class')
        }
        else{
            M1 = glm(yy~., data = data.frame(yy=y, x[,feature]), family = binomial(link="logit"))
            yhat = levels(y)[(predict(M1, type = 'response')>0.5)+1]
        }
        #tree$value = M1
        left_index = (yhat == levels(y)[1])
        if(sum(left_index) == n || sum(left_index) == 0){
            tree$value = levels(y)[1+(mean(y == levels(y)[2])>0.5)]
            return(tree)
        }else{
            beta = coef(M1)
            beta[is.na(beta)] = 0
            tree$value = beta
            tree$left = Tree_build(x = x[left_index,], y = y[left_index], mtry = mtry, nodesize = nodesize, method = method, ridge = ridge)
            tree$right = Tree_build(x = x[!left_index,], y = y[!left_index], mtry = mtry, nodesize = nodesize, method = method, ridge = ridge)
            return(tree)
        }
    }
}


#logit function prediction
logit_predict = function(beta, x, level){
    if(is.null(dim(x))){
        n = 1
    }else{
        n = dim(x)
    }

    link = x%*%beta[-1] + beta[1]
    y = exp(link)/(1+exp(link))
    y[link>50] = 1
    y[link< -50] = 0
    #if(any(is.na(y)))
    #    print(c(link, beta, x))
    return(level[1+(y>0.5)])
}

#prediction function for a single logistic tree
predict_tree = function(tree, newdata, level){
    if(is.null(dim(newdata))){
        n = 1
    }else{
        n = dim(newdata)[1]
    }
    if((is.null(tree$left))||(is.null(tree$right))){
        predicted = matrix(tree$value, n, 1)
        return(predicted)
    }else{

        predicted = matrix(0, n, 1)
        if(n == 1){
            direct = logit_predict(tree$value, x = t(as.matrix(newdata[tree$feature])), level = level)
        }else{
            direct = logit_predict(tree$value, x = newdata[,tree$feature], level = level)
        }

        if(n == 1 && direct == level[1])
            predicted[1] = predict_tree(tree$left, newdata, level = level)
        else if(n == 1)
            predicted[1] = predict_tree(tree$right, newdata, level = level)
        else{       
            direct_left = (direct == level[1])
            if(sum(direct_left)>0)
                predicted[direct_left] = predict_tree(tree$left, newdata[direct_left,], level = level)
            if(sum(direct_left)<n)
                predicted[!direct_left] = predict_tree(tree$right, newdata[!direct_left,], level = level)
        }
        return(predicted)
    }
}

