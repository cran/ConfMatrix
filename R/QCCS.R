#' @title Quality Control Columns Set
#' @description The difference between a QCCS and a confusion matrix is
#' that while forming a confusion matrix requires that the reference and
#' the product be more or less equivalent, for the QCCS it is required
#' that the reference be actually of higher quality than the product.
#' This forces us to leave the marginals corresponding to the reference
#' fixed. That is why we work by columns. In this way, the QCCS class
#' works with a confusion matrix expressed as a set of column vectors
#' and it will be analyzed by columns. A QCCS is constructed by comparing
#' a sample of a set of common positions in the product and the ground
#' truth. Appropriate sampling methods must be applied to generate the
#' QCCS. It is considered that the classes of the ground truth
#' correspond to the columns and that the classes of the
#' product to be valued correspond to the rows. On the other hand, the
#' concept of QCCS is directly linked to quality control, so the
#' specifications of this control must be indicated \insertCite{QCCS}{ConfMatrix}.
#' Specifications are stated as percentages. E.g. for class "A" under
#' consideration, a minimum quality value is established (e.g. better than 90%), and
#' maximum values of confusion with other categories (e.g. confusion
#' between A and B less than 5%). The specifications are proportions of
#' a multinomial. First, an object of this class of object must be
#' created (instantiated) and then the methods that offer the index
#' calculations will be invoked.
#' @export QCCS
#' @note  Error Messages: List of possible errors:
#' \itemize{
#'  \item \code{Error type 1}: Different number of data vectors and probability.
#'  \item \code{Error type 2}: Different number of elements in the pair of data
#'  vectors and probabilities.
#'  \item \code{Error type 3}: The sum of the elements of the data vectors is 0.
#'  \item \code{Error type 4}: The sum of each probability vectors must be 1.
#'  \item \code{Error type 5}: Some element of the data vector is negative.
#'  \item \code{Error type 6}: Some element of the probability vector is negative.
#'}
#' @references
#' \insertRef{alba2020}{ConfMatrix}
#'
#' \insertRef{QCCS}{ConfMatrix}
#' @importFrom R6 R6Class
#' @importFrom stats dmultinom pchisq
#' @importFrom Rdpack reprompt
#'
#'
#' @aliases QCCS


QCCS <- R6Class("QCCS",
  cloneable=FALSE,
   public = list(
    #' @field Vectors
    #'\verb{
    #'List of integer values data for the vectors.
    #'}
    Vectors = NULL,
    #' @field Prob
    #'\verb{
    #'List of probability values corresponding to each of the vectors.
    #'}
    Prob = NULL,
    #' @field ID
    #'\verb{
    #'Identifier. It is a character string with a maximum length of 50 characters.
    #'By default,} \eqn{QCCS_i} \verb{will be taken as identification. Where} \eqn{i \in [1,999]} \verb{will be the
    #'number of QCCS instances already defined in the session.
    #'}
    ID=NULL,
    #' @field Date
    #'\verb{
    #'Date provided by the user in format DDMMYYYY, "DD-MM-YYYY", "DD/MM/YYYY".
    #'By default the date provided by the system will be taken.
    #'}
    Date=NULL,
    #' @field ClassNames
    #'\verb{
    #' Name of the classes. It is given by a character strings vector whose elements
    #' are the name of the classes. Each element of the vector is a string of maximum
    #' 20 characters. By default for the column elements they will be} \eqn{PC_i'} \verb{ (Producer
    #' class).}
    #'
    ClassNames=NULL,
    #' @field Source
    #'\verb{
    #' Indicates where the "vectors" and "prob" parameters come from (article, project,
    #' etc.). It is suggested to enter a reference or a DOI. A character string with
    #' a maximum length of 80 characters can be entered. By default, is NULL.
    #'}
    Source=NULL,


    #' @description Public method to create an instance of the QCCS class.
    #' At the time of creation, column set data and specification values
    #' must be provided. The same number of data and as specification values
    #' must be entered, and the pairs of data-specifications vectors must
    #' have the same size, otherwise an error will be provided.
    #' The optional possibility of adding metadata to the matrix is offered.
    #' The values of the data vectors represent the classes of ground truth.
    #' @param Vectors
    #' \verb{
    #' List of integer values data for the vectors.
    #' }
    #' @param Prob
    #' \verb{
    #' List of probability values corresponding to each of the vectors.
    #' }
    #' @param ID
    #'\verb{
    #'Identifier. It is a character string with a maximum length of 50 characters.
    #'By default,} \eqn{QCCS_i} \verb{will be taken as identification. Where} \eqn{i \in [1,999]} \verb{will be
    #'the number of QCCS instances already defined in the session.
    #'}
    #' @param Date
    #'\verb{
    #' Date provided by the user in format DDMMYYYY, "DD-MM-YYYY", "DD/MM/YYYY".
    #' By default the date provided by the system will be taken.
    #'}
    #' @param ClassNames
    #' \verb{
    #' Name of the classes. It is given by a character strings vector whose elements
    #' are the name of the classes. Each element of the vector is a string of maximum
    #' 20 characters. By default for the column elements they will be} \eqn{PC_i'} \verb{ (Producer
    #' class).}
    #'
    #' @param Source
    #' \verb{
    #' Indicates where the "vectors" and "prob" parameters come from (article, proj-
    #' ect, etc.). It is suggested to enter a reference or a DOI. A character string
    #' with a maximum length of 80 characters can be entered. By default, is NULL.
    #' }
    #' @examples
    #' Vectors<-list(c(47,4,0),c(44,5,3))
    #' Prob<-list(c(0.95,0.04,0.01),c(0.88,0.1,0.02))
    #' A<-QCCS$new(Vectors,Prob,
    #' Source="Ariza-Lopez et al. 2019")
    #'
    #' @aliases NULL

  initialize = function(Vectors,Prob,ID=NULL,Date=NULL,ClassNames=NULL,Source=NULL) {


# Optional values ---------------------------------------------------------


    self$Vectors <- Vectors
    self$Prob <- Prob

    if(is.null(ID)){
      sequence<- private$sequence()
      self$ID <- paste("QCCS_",sequence,sep="")
    }else{
      self$ID<-substr(ID,1,50)
    }
    if(!is.null(Date)){
      self$Date<-Date
    }else{self$Date <- Sys.Date()}

    colname<-c()
    if (!is.null(ClassNames)) {
      self$ClassNames <- ClassNames
      for (i in 1:length(self$Vectors)) {
        colname <- c(colname, sprintf("PC_%.20s", self$ClassNames[i]))
      }
      names(self$Vectors) <- colname

    } else {
      self$ClassNames <- ClassNames
      for (i in 1:length(self$Vectors)) {
        colname <- c(colname, sprintf("PC_%d", i))
      }
      names(self$Vectors) <- colname
    }
    if(!is.null(Source)){
      self$Source <- substr(Source,1,80)
    }else{self$Source<-NULL}

# Initializing values and checking if they are correct values -----------

    n <- length(self$Vectors)
    m <- length(self$Prob)
      if (n != m) {
        stop("Error type 1: There must be the same number of\ndata columns as Probability columns")
      }

      for (i in 1:n) {
        vi <- self$Vectors[[i]]
        pi <- self$Prob[[i]]
        ni <- length(vi)
        mi <- length(pi)

        if ((ni != mi) == TRUE) {
          stop("Error type 2: The vectors and their corresponding\nprobabilities must have the same size\n")
        }

        if(sum(vi)==0){
          stop("Error type 3: The sum of the elements of the\ndata vectors is 0\n")
        }

        if(sum(pi)!=1){
          stop("Error type 4: The sum of each probability\nvectors must be 1\n")
        }

        if(length(vi[vi<0])>0){
          stop("Error type 5: Some element of the data\nvector is negative\n")
        }

        if(length(pi[pi<0])>0){
          stop("Error type 6: Some element of the probability\nvector is negative.\n")
        }
      }
  },






# print function ----------------------------------------------------------

      #' @description Public method that shows all the data entered
      #' by the user.
      #' @return QCCS object identifier, Date, name of classes, source
      #' of data and data vectors and probability.
      #' @examples
      #' Vectors<-list(c(18,0,3,0),c(27,19))
      #' Prob<-list(c(0.85,0.1,0.03,0.02),c(0.8,0.2))
      #' A<-QCCS$new(Vectors,Prob,
      #' Source="Alba-Fernández et al. 2020")
      #' A$print()
      #'
      #' @aliases NULL

    print=function(){
      cat("Identifier (ID)\n", self$ID, "\n")
      cat("-------------------------------------\n")
      cat(sprintf("Date\n %s \n", self$Date))
      cat("-------------------------------------\n")
      cat("Source\n", self$Source, "\n")
      cat("-------------------------------------\n")
      for(i in 1:length(self$Vectors)){
        cat("Name of Class|",names(self$Vectors)[i], "\n")
        cat("Vector       |",self$Vectors[[i]],"\n")
        cat("Probability  |",self$Prob[[i]],"\n")
        cat("-------------------------------------\n")

      }

    },

# Multinomial or binomial Exact Tests -------------------------------------


      #' @description Public method that using a QCCS object
      #' instance calculates whether the data meets specifications.
      #' An exact test is applied to each of the multinomials
      #' that are defined for each column.
      #' The Bonferroni method is used.
      #' The references \insertCite{QCCS}{ConfMatrix} and \insertCite{alba2020}{ConfMatrix}
      #' are followed for the computations.
      #' @param a \verb{
      #' significance level. By default a=0.05.
      #' }
      #' @return A list of the "htest" class containing the results of the hypothesis test.
      #' The p-value returned is the lowest of those obtained for the data analyzed.
      #' In addition, the Bonferroni criterion value, the p-values obtained for each column,
      #' the original data vectors and the probability vectors are also returned as
      #' parameters of the htest class.
      #' @examples
      #' \donttest{
      #' Vectors<-list(c(47,4,0),c(40,5,3))
      #' Prob<-list(c(0.95,0.04,0.01),c(0.88,0.1,0.02))
      #' A<-QCCS$new(Vectors,Prob,
      #' Source="Ariza-Lopez et al. 2019")
      #' A$Exact.test()
      #' }
      #'
      #' @aliases NULL


    Exact.test = function(a=NULL) {
      if(is.null(a)){
        a<-0.05
      }else{a<-a}

      n <- length(self$Vectors)
      m <- length(self$Prob)
      sol <- c()
      p_value<-c()
      p_value1<-0
        for (i in 1:n) {
          vi <- self$Vectors[[i]]
          pi <- self$Prob[[i]]
          ni <- length(vi)
          mi <- length(pi)
          s <- sum(vi)

          if(ni==2){
            p_value1<-private$test.2tol(vi,pi)[1]
          } else
          if(ni!=2){
            p_value1<-private$test.ntol(vi,pi)$p.valor
          }
          p_value<-c(p_value,p_value1)
        }

     sol <- c(sol, p_value)

      htest_result <- list(
        method = "Exact Test with Bonferroni Correction",
        p.value = min(p_value),
        data.name = "Vectors and Probabilities",
        Bonferroni.criterion=a/length(p_value),
        Individual.pvalues=p_value,
        OriginalVectors=self$Vectors,
        OriginalProb=self$Prob
      )

      class(htest_result) <- "htest"
      return(htest_result)
    },


# ji multinomial or binomial test ----------------------------------------------

      #' @description Public method that using a QCCS object instance
      #' calculates whether the data meets specifications in each of the classes.
      #' The Chi square test is used. The Bonferroni method is used.
      #' The references \insertCite{QCCS}{ConfMatrix} and \insertCite{alba2020}{ConfMatrix}
      #' are followed for the computations.
      #' @param a \verb{
      #' significance level. By default a=0.05.
      #' }
      #' @return A list of the "htest" class containing the results of the hypothesis test.
      #' The p-value returned is the lowest of those obtained for the data analyzed.
      #' In addition, the Bonferroni criterion value, the obtained p-values, the degrees of
      #' freedom and the statistics obtained for each column, the original data vectors
      #' and the probability vectors are also returned as parameters of the htest class.
      #' @examples
      #' Vectors<-list(c(18,0,3,0),c(27,19))
      #' Prob<-list(c(0.85,0.1,0.03,0.02),c(0.8,0.2))
      #' A <- QCCS$new(Vectors,Prob,
      #' Source="Alba-Fernández et al. 2020")
      #' A$Ji.test()
      #'
      #' @aliases NULL

    Ji.test=function(a=NULL){
      if(is.null(a)){
      a<-0.05
      }else{a<-a}

      n <- length(self$Vectors)
      m <- length(self$Prob)
      p_value<-c()
      k<-0
      statistics <- c()
      df <- c()
        for (j in 1:n) {
        vi <- self$Vectors[[j]]
        pi <- self$Prob[[j]]
        ni <- length(vi)
        mi <- length(pi)
        k<-ni-1
        pj<-c()
        for (i in 1:ni){
          pj<-c(pj,(vi[i]-sum(vi)*pi[i])/sqrt(sum(vi)*pi[i]))
        }
        ST<-sum(pj^2)
        statistics <- c(statistics, ST)
        df <- c(df, k)
        pvalue<-pchisq(ST, k, lower.tail=FALSE)
        p_value<-c(p_value,pvalue)
        }

        htest_result <- list(
              p.value = min(p_value),
              method = "Chi-squared test with Bonferroni method",
              data.name = "Vectors and Probabilities",
              Bonferroni.criterion=a/length(p_value),
              Individual.pvalues=p_value,
              Individual.statistics = statistics,
              Individual.parameters = df,
              OriginalVectors=self$Vectors,
              OriginalProb=self$Prob
        )

        class(htest_result) <- "htest"
        return(htest_result)
    },


# ji global multinomial or binomial test ----------------------------------------------

      #' @description Public method that using a QCCS
      #' object instance calculates whether the data meets
      #' specifications. The Chi square test is used.
      #' The references
      #' \insertCite{QCCS}{ConfMatrix} and \insertCite{alba2020}{ConfMatrix}
      #' are followed for the computations.
      #' @param a \verb{
      #' significance level. By default a=0.05.
      #' }
      #' @return A list of class "htest" containing the results of the hypothesis test.
      #' In addition, the original data vectors and the probability vectors are also
      #' returned.
      #' @examples
      #' Vectors<-list(c(18,0,3,0),c(27,19))
      #' Prob<-list(c(0.85,0.1,0.03,0.02),c(0.8,0.2))
      #' A <- QCCS$new(Vectors,Prob,
      #' Source="Alba-Fernández et al. 2020")
      #' A$JiGlobal.test()
      #'
      #' @aliases NULL

    JiGlobal.test=function(a=NULL){
      if(is.null(a)){
        a<-0.05
      }else{a<-a}

    n <- length(self$Vectors)
    m <- length(self$Prob)
    p_value<-c()
    Suma<-0
    S<-list()
    k<-0
      for (j in 1:n) {
      vi <- self$Vectors[[j]]
      pi <- self$Prob[[j]]
      ni <- length(vi)
      mi <- length(pi)
      k<-k+(ni-1)
      pj<-c()
        for (i in 1:ni){
          pj<-c(pj,(vi[i]-sum(vi)*pi[i])/sqrt(sum(vi)*pi[i]))
        }
      S[[j]]<-pj

      ST<-pj^2

      Suma<-Suma+sum(ST)
      }


    p_value<-pchisq(Suma, k, lower.tail=FALSE)

    htest_result <- list(
        statistic = c(X2 = Suma),
        parameter = c(df = k),
        p.value = p_value,
        method = "Global Chi-squared test",
        data.name = "Vectors and Probabilities",
        OriginalVectors = self$Vectors,
        OriginalProb = self$Prob
      )

    class(htest_result) <- "htest"
    return(htest_result)
    }



  ),



# Private functions -------------------------------------------------------



  private = list(

    combines = function(n, N) {
      result <- NULL
      stack <- list(list(comb = NULL, remaining = N, index = 1))
      pb <- txtProgressBar(min = 0,
                           max = choose(N+1+n-1,n),
                           style = 3,
                           width = 50,
                           char = "=")
      k<-0
      while (length(stack) > 0) {
        current <- stack[[length(stack)]]
        stack <- stack[-length(stack)]

        comb <- current$comb
        remaining <- current$remaining
        index <- current$index
          if (length(comb) == n) {
            k=k+1
            setTxtProgressBar(pb, k)
            result <- rbind(result, comb)
          } else {
            for (i in 0:remaining) {
            stack <- c(stack, list(list(comb = c(comb, i),
                     remaining = remaining - i, index = index + 1)))
          }
        }
      }
      close(pb)
      return(result)
    },

    sequence = function() {
      variables <- ls(envir = .GlobalEnv)
      es_qccs <- sapply(variables, function(x) {
        obj <- tryCatch(get(x, envir = .GlobalEnv), error = function(e) NULL)
        if (!is.null(obj)) {
          return(inherits(obj, "QCCS"))
        } else {
          return(FALSE)
        }
      }, USE.NAMES = FALSE)

      number_qccs <- sum(es_qccs, na.rm = TRUE)
      return(number_qccs + 1)
    },

# Functions depending on the number of tolerances -------------------------

    test.ntol = function(M, p){
      n<-length(M)
      N <- sum(M)

      Samples1<-NULL
      Samples<-private$combines(n,N)
      Samples1<- matrix(ncol = n,nrow=0)

        for (i in 1:dim(Samples)[1]) {
          if (sum(Samples[i,]) == N) {
          Samples1<-rbind(Samples1,Samples[i,])
          }
        }
      Q <- matrix(ncol = n+2, nrow = dim(Samples1)[1])

      for (i in 1:dim(Samples1)[1]) {
        Q[,1:n]<-Samples1[,1:n]
        Q[i,n+1]<-dmultinom(Samples1[i,1:n],prob=p)
        Q[1,n+2]<-Q[1,n+1]
        if(i>1){
          Q[i,n+2]<-Q[i, n+1] + Q[(i - 1), n+2]
        }
      }
      A <- matrix(nrow = 0, ncol = n+2)
      cond<-"Q[j,1]==M[1]"
      cond1<-paste(" Q[j,", 2:n, "] == M[", 2:n, "]", sep = "")
      cond<-c(cond,cond1)
      condition<-noquote(cond)

      for(j in 1:dim(Q)[1]) {
        if ( (Q[j, 1] < M[1])){
          A <- rbind(A, Q[j, ])
        }

        for(ni in 1:n-1){
        condition1<-paste(cond[1:ni], collapse = " & ")

          if(( eval(parse(text=condition1)))){
            if(Q[j, ni+1] < M[ni+1]){
              A<-rbind(A,Q[j,])
            }
          }
        }
        if(Q[j,1]==M[1] & Q[j,2]==M[2] & Q[j,3]==M[3] ){
              A<-rbind(A,Q[j,])
        }
      }

      p_valor<-sum(A[,n+1])
      results <- list(A = A, p.valor = p_valor)

    return(results)
    },


# Function for 2 tolerances -----------------------------------------------



    test.2tol=function(M, p){
      N<-sum(M)
      p<-pbinom(M,N,p)
      p_value<-p[1]
     return(p_value)
    }

 )
)

