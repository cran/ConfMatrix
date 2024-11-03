#' @title Confusion matrix
#' @description
#' The ConfMatrix class works with confusion matrices, thus providing
#' the possibility of calculating several indices with their
#' corresponding variances and confidence intervals. A confusion matrix
#' is constructed by comparing a sample of a set of common positions in
#' the product and the ground truth. Appropriate sampling methods must
#' be applied to generate the confusion matrix. It is considered that
#' the classes of the ground truth correspond to the columns
#' and that the classes of the product to be valued correspond
#' to the rows. First, an object of this class of object must be created
#' (instantiated) and then the methods that offer the index calculations
#' will be invoked. Mnemonic method names are proposed and are therefore
#' long, for example methods that provide averages start with "AV" and
#' those that provide combinations start with "Comb". Methods related
#' to a specific thematics class end with the ending "_i".
#'
#' @note Error Messages: List of possible errors:
#' \itemize{
#'  \item \code{Error type 1}: Non-square matrix.
#'  \item \code{Error type 2}: Single element matrix.
#'  \item \code{Error type 3}: Negative values.
#'  \item \code{Error type 4}: Sum of elements 0.
#'  \item \code{Error type 5}: Sum of rows 0.
#'  \item \code{Error type 6}: Sum of columns 0.
#'  \item \code{Error type 7}: It is not a matrix.
#'}
#' @export ConfMatrix
#' @importFrom R6 R6Class
#' @importFrom Rdpack reprompt
#' @importFrom gridExtra grid.arrange
#' @importFrom ggplot2 ggplot aes xlab ylab theme geom_point geom_errorbar
#' @references
#'
#' \insertRef{alba2020}{ConfMatrix}
#'
#' \insertRef{ariza2011}{ConfMatrix}
#'
#' \insertRef{book}{ConfMatrix}
#'
#' \insertRef{cohen1960}{ConfMatrix}
#'
#' \insertRef{congalton2008}{ConfMatrix}
#'
#' \insertRef{condkappa}{ConfMatrix}
#'
#' \insertRef{fienberg1970}{ConfMatrix}
#'
#' \insertRef{finn1993}{ConfMatrix}
#'
#' \insertRef{fleiss1969}{ConfMatrix}
#'
#' \insertRef{foody1992}{ConfMatrix}
#'
#' \insertRef{garcia2018}{ConfMatrix}
#'
#' \insertRef{ghosh2002}{ConfMatrix}
#'
#' \insertRef{goodman1968analysis}{ConfMatrix}
#'
#' \insertRef{hellden1980}{ConfMatrix}
#'
#' \insertRef{koukoulas2001}{ConfMatrix}
#'
#' \insertRef{labatut2011}{ConfMatrix}
#'
#' \insertRef{liu2007}{ConfMatrix}
#'
#' \insertRef{ma1995Tau}{ConfMatrix}
#'
#' \insertRef{munoz2016}{ConfMatrix}
#'
#' \insertRef{naesset1996}{ConfMatrix}
#'
#' \insertRef{pontius2014}{ConfMatrix}
#'
#' \insertRef{diffeR}{ConfMatrix}
#'
#' \insertRef{rosenfield1986}{ConfMatrix}
#'
#' \insertRef{short1982}{ConfMatrix}
#'
#' \insertRef{strehl2002relationship}{ConfMatrix}
#'
#' \insertRef{strehl2002}{ConfMatrix}
#'
#' \insertRef{tung1988}{ConfMatrix}
#'
#' \insertRef{turk1979gt}{ConfMatrix}
#'
#' \insertRef{turk2002}{ConfMatrix}
#'
#'
#' @section Mathematical elements:
#' \itemize{
#'   \item \eqn{x_{ii}}: diagonal element of the matrix.
#'   \item \eqn{x_{ij}}: element \eqn{i,j} of the matrix.
#'   \item \eqn{x_{i+}}: sum of all elements in rows \eqn{i}.
#'   \item \eqn{x_{+j}}: sum of all elements in column \eqn{j}.
#'   \item \eqn{M}: number of classes.
#'   \item \eqn{\overline{x}_{i+}}: sum of all elements of row \eqn{i} except element \eqn{i} of the diagonal.
#'   \item \eqn{\overline{x}_{+i}}: sum of all elements of column \eqn{i} except element \eqn{i} of the diagonal.
#'   \item \eqn{N_{Total}}:Total count of elements in the instance's Confusion Matrix.
#'   \deqn{N_{Total}=\sum_{i,j}^M x_{ij}}
#'   \item \eqn{N_i / N_j}: Total count of elements in row \eqn{i} or column \eqn{j}.
#'   \deqn{N_{i}=x_{i+}}
#'   \deqn{N_{j}=x_{+j}}
#'   \item \eqn{N_{ij}}: Total count of elements in row \eqn{i} and column \eqn{j}.
#'   \deqn{N_{ij}=x_{i+}+x_{+j}-x_{ii}}
#'}
#'
#'
#' @aliases ConfMatrix




ConfMatrix <- R6Class("ConfMatrix",
  cloneable=FALSE,
   public = list(
    #' @field Values
    #'\verb{
    #'Matrix of integer values. An matrix must be added.
    #'}
    Values = NULL,
    #' @field ID
    #'\verb{
    #' Identifier. It is a character string with a maximum length of 50 characters.
    #' By default,} \eqn{CM_i} \verb{will be taken as identification. Where} \eqn{i \in [1,999]} \verb{will be the
    #' number of ConfMatrix instances already defined in the session.
    #' }
    ID = NULL,
    #' @field Date
    #' \verb{
    #' Date provided by the user in format DDMMYYYY, "DD-MM-YYYY", "DD/MM/YYYY".
    #' By default the date provided by the system will be taken.
    #'
    #'}
    Date = NULL,
    #' @field ClassNames
    #'\verb{
    #' Name of the classes. It is given by a character strings vector whose elements
    #' are the name of the classes. Each element of the vector is a string of maximum
    #' 20 characters. By default for the column elements they will be} \eqn{PC_i}
    #' (Producer class) \verb{and for the elements of row} \eqn{UC_i}\verb{ (User class), with} \eqn{i} \verb{being the correspond-
    #' ing row or column number.
    #' }
    ClassNames=NULL,
    #' @field Source
    #'  \verb{
    #' Indicates where the matrix comes from (article, project, etc.). It is suggest-
    #' ed to enter a reference or a DOI. A character string with a maximum length of
    #' 80 characters can be entered. By default, is NULL.
    #' }
    Source=NULL,


    #' @description
    #' Public method to create an instance of the ConfMatrix class.
    #' When creating it, values must be given to the matrix. The values
    #' of the matrix must be organized in such a way that the columns
    #' represent the classes in the reference and the rows represent
    #' the classes in the product being evaluated. The creation of a
    #' ConfMatrix instance includes a series of checks on the data. If
    #' checks are not met, the system generates coded error messages.
    #' The optional possibility of adding metadata to the matrix is offered.
    #' @param Values
    #'\verb{
    #' Matrix of integer values. A matrix must be added.
    #'}
    #'
    #' @param ID
    #'\verb{
    #' Identifier. It is a character string with a maximum length of 50 characters.
    #' By default,} \eqn{CM_i} \verb{will be taken as identification. Where} \eqn{i \in [1,999]} \verb{will be the
    #' number of ConfMatrix instances already defined in the session.
    #'}
    #' @param Date
    #'\verb{
    #' Date provided by the user in format DDMMYYYY, "DD-MM-YYYY", "DD/MM/YYYY".
    #' By default the date provided by the system will be taken.
    #'
    #'}
    #' @param ClassNames
    #' \verb{
    #' Name of the classes. It is given by a character strings vector whose elements
    #' are the name of the classes. Each element of the vector is a string of maximum
    #' 20 characters. By default for the column elements they will be} \eqn{PC_i}
    #' (Producer class) \verb{and for the elements of row} \eqn{UC_i}\verb{ (User class), with} \eqn{i} \verb{being the correspond-
    #' ing row or column number.
    #' }
    #' @param Source
    #' \verb{
    #' Indicates where the matrix comes from (article, project, etc.). It is suggest-
    #' ed to enter a reference or a DOI. A character string with a maximum length of
    #' 80 characters can be entered. By default, is NULL.
    #' }
    #'
    #' @return Object of the ConfMatrix class, or an error message.
    #'
    #' @examples
    #' A<-matrix(c(65,6,0,4,4,81,11,7,22,5,85,3,24,8,19,90),nrow=4,ncol=4)
    #' cm<-ConfMatrix$new (A,ID="5",Date="27-10-2023",Source="Congalton and Green,
    #' 2008")
    #'
    #' @aliases NULL

initialize = function(Values,ID=NULL,Date=NULL,ClassNames=NULL,Source=NULL) {


# Initializing values -----------------------------------------------------


  self$Values<-Values

  if(is.null(ID)){
    sequence <- private$sequence()
    self$ID <- paste("CM_", sequence, sep="")
  }else{
    self$ID<-substr(ID,1,50)
  }

  if(!is.null(Date)){
    self$Date<-Date
  }else{self$Date <- Sys.Date()}

  colname<-c()
  rowname<-c()
  if(!is.null(ClassNames)){
    self$ClassNames<-ClassNames

    for(i in 1:sqrt(length(self$Values))){
      colname<-c(colname,sprintf("PC_%.20s",self$ClassNames[i]))
      rowname<-c(rowname,sprintf("UC_%.20s",self$ClassNames[i]))
    }

    self$ClassNames <- ClassNames
    colnames(self$Values)<-colname
    rownames(self$Values)<-rowname

  }else{
    for(i in 1:nrow(self$Values)){
      colname<-c(colname,sprintf("PC_%d",i))
      rowname<-c(rowname,sprintf("UC_%d",i))
    }
    colnames(self$Values)<-colname
    rownames(self$Values)<-rowname
    self$ClassNames<-c(1:nrow(self$Values))
      }


  if(!is.null(Source)){
    self$Source <- substr(Source,1,80)
  }else{self$Source<-NULL}

  nk<-nrow(self$Values)
  nfilas <- nrow(self$Values)
  ncolumnas <- ncol(self$Values)

# Matrix check ------------------------------------------------------------

   if(nfilas != ncolumnas){
     stop("Error type 1: Non-square matrix\n")
   }

   if(nk==1){
     stop("Error type 2: Single element matrix\n")
   }

  if(length(self$Values[self$Values<0])>0){
    stop("Error type 3: negative values.")
  }

  if(sum(self$Values)==0 ){
     stop("Error type 4: Sum of elements 0\n")
  }
  
   if(sum(apply(self$Values,1,sum))==0 ){
     stop("Error type 5: Sum of rows 0\n")
   }
  
   if(sum(apply(self$Values,2,sum))==0 ){
     stop("Error type 6: Sum of columns 0\n")
   }
  
   if(is.matrix(self$Values) == FALSE){
     stop("Error type 7: It is not a matrix\n")
   }
  
  },

# graphics ----------------------------------------------------------------


      #' @description Public method that provides a graph of the indices of
      #' the functions ConfMatrix$OverallAcc, ConfMatrix$Kappa,
      #' ConfMatrix$Tau, ConfMatrix$AvHellAcc and ConfMatrix$AvShortAcc
      #'  with their corresponding standard deviations.
      #' @return A graph of the indices of the functions OverallAcc, Kappa,
      #' Tau, AvHellAcc, AvShortAcc with their corresponding
      #' standard deviations.
      #' @examples
      #' A<-matrix(c(65,6,0,4,4,81,11,7,22,5,85,3,24,8,19,90), nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A,Source="Congalton and Green 2008")
      #' p$plot.index()
      #'
      #' @aliases NULL

    plot.index=function(){
      index1<-c(self$OverallAcc()[[1]],self$Tau()[[1]],self$Kappa()[[1]],
      self$AvHellAcc()[[1]],self$AvShortAcc()[[1]])

      var<-c(self$OverallAcc()[[2]],self$Tau()[[2]],self$Kappa()[[2]],
             self$AvHellAcc()[[2]],self$AvShortAcc()[[2]])

      desv<-sqrt(var)

      index<-c("OverallAcc","Tau","Kappa","AvHellAcc","AvShortAcc")
      datos<-data.frame(index1,desv,index)

      graf<-ggplot(datos,aes(x=index,y=index1,colour=index))
      graf<-graf+geom_point(size=3)+ylab("values")+xlab("Global Indices")
      graf<-graf+geom_errorbar(aes(ymin=index1-desv, ymax=index1+desv),
            width = 0.5)
    return(graf)
    },



      #' @description Public method that provides a graph for the user’s
      #' and producer’s accuracies and standard deviations.
      #' @return The graph of the accuracy index of users and producers
      #' with their corresponding standard desviation.
      #' @examples
      #' A<-matrix(c(65,6,0,4,4,81,11,7,22,5,85,3,24,8,19,90), nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A,Source="Congalton and Green 2008")
      #' p$plot.UserProdAcc()
      #'
      #' @aliases NULL

    plot.UserProdAcc=function(){

      index1<-c(self$UserAcc()[[1]])
      index2<-c(self$ProdAcc()[[1]])
      var<-c(self$UserAcc()[[2]])
      var2<-c(self$ProdAcc()[[2]])
      desv<-sqrt(var)
      desv2<-sqrt(var2)

      index<-c()
      for (i in 1:length(self$UserAcc()[[1]])) {
        index<-c(index,sprintf("Class %d",i))
      }

      datos<-data.frame(index1,desv,index)
      datos1<-data.frame(index2,desv2,index)
      graf<-ggplot(datos,aes(x=index,y=index1,colour=index))
      graf<-graf+geom_point(size=3)+ylab("values")+xlab("User Accuracy")
      graf<-graf+geom_errorbar(aes(ymin=index1-desv, ymax=index1+desv), width = 0.5)+theme(legend.position = "none")

      graf1<-ggplot(datos1,aes(x=index,y=index2,colour=index))
      graf1<-graf1+geom_point(size=3)+ylab("values")+xlab("Producer Accuracy")
      graf1<-graf1+geom_errorbar(aes(ymin=index2-desv2, ymax=index2+desv2), width = 0.5)+theme(legend.position = "none")

    return(grid.arrange(graf, graf1, ncol = 2))
    },

# print function ----------------------------------------------------------

      #' @description Public method that shows all the data entered
      #' by the user for a instance.
      #' @return ConfMatrix object identifier, date, class name, data
      #' source and confusion matrix.
      #' @examples
      #' A<-matrix(c(65,6,0,4,4,81,11,7,22,5,85,3,24,8,19,90), nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A,ClassNames=c("Deciduous","conifer","agriculture",
      #' "shrub"),Source="Congalton and Green 2008")
      #' p$print()
      #'
      #' @aliases NULL

    print=function(){
      cat("Identifier (ID)\n", self$ID, "\n")
      cat("-------------------------------------\n")
      cat(sprintf("Date\n %s \n", self$Date))
      cat("-------------------------------------\n")
      cat("ClassNames\n", self$ClassNames, "\n")
      cat("-------------------------------------\n")
      cat("Source\n", self$Source, "\n")
      cat("-------------------------------------\n")
      cat("Confusion Matrix\n")
      print(self$Values)
    },


# All parameters ----------------------------------------------------------


      #' @description Public method in which multiple parameters are
      #' calculated for the given confusion matrix. This method is
      #' equivalent to ConfMatrix$OverallAcc,ConfMatrix$UserAcc,
      #' ConfMatrix$ProdAcc,ConfMatrix$Kappa and ConfMatrix$MPseudoZeroes.
      #' @return The following list of elements: the confusion matrix,
      #' dimension, total sum of cell values, overall accuracy, overall
      #' accuracy variance, global kappa index, global kappa simplified
      #' variance, producer accuracy by class, user accuracy by class,
      #' and pseudoceros matrix.
      #' @examples
      #' A<-matrix(c(65,6,0,4,4,81,11,7,22,5,85,3,24,8,19,90), nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A,Source="Congalton and Green 2008")
      #' p$AllParameters()
      #'
      #' @aliases NULL

     AllParameters=function(){
      OverallAcc<-self$OverallAcc()
      dimension <- nrow(self$Values)
      SumaMatriz <-sum(self$Values)
      PAcuerdo <- OverallAcc[[1]]
      ExProdu <- self$ProdAcc()[[1]]
      UserAcc <-self$UserAcc()[[1]]
      PAAzar <- sum((private$sumfil(self$Values)*private$sumcol(self$Values)))/(SumaMatriz*SumaMatriz)
      Kappa <- self$Kappa()[[1]]
      VarPAcuerdo <- self$OverallAcc()[[2]]
      VarKappa <-  self$Kappa()[[2]]
      MPseudoceros <- self$MPseudoZeroes()[[2]]
      salida<-list(Matrix=self$Values, Dimension =dimension, n=SumaMatriz,
                   OverallAcc=PAcuerdo, VarOverallAcc=VarPAcuerdo, Kappa=Kappa,
                   VarKappa=VarKappa,ProdAcc=ExProdu,UserAcc=UserAcc,
                   MPseudoceros=MPseudoceros)
     return(salida)
     },


# UserAcc -----------------------------------------------------------------



      #' @description Public method for deriving the index called user’s
      #' accuracy for all the classes in a ConfMatrix object instance.
      #' The user's accuracy for the class \eqn{i} of a thematic product is
      #' calculated by dividing the value in the diagonal of class \eqn{i} by
      #' the sum of all values in the row of the class \eqn{i} (row marginal).
      #' The method also offers the variance and confidence interval.
      #' The reference \insertCite{congalton2008;textual}{ConfMatrix} is followed
      #' for the computations.
      #' \deqn{
      #' UserAcc=\dfrac{x_{ii}}{x_{i+}}
      #' }
      #'  \deqn{
      #' \sigma^2_{UserAcc}=\dfrac{UserAcc \cdot (1-UserAcc)}{N_{i}}
      #' }
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #' @return A list of vectors, containing the user’s accuracy real values for
      #' all classes, their variances and confidence intervals for each class,
      #' respectively.
      #'
      #' @examples
      #' A<-matrix(c(65,6,0,4,4,81,11,7,22,5,85,3,24,8,19,90), nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A,Source="Congalton and Green 2008")
      #' p$UserAcc()
      #'
      #' @aliases NULL



     UserAcc = function(a=NULL){
     n <- sqrt(length(self$Values))
     UserAcc <- rep(0,n)
     VarUserAcc<-rep(0,n)
     ConfInt<-list()
       for (i in 1:n){
         UserAcc[i] <- self$Values[i,i] / private$sumfil(self$Values)[i]
         VarUserAcc[i]<-abs((UserAcc[i]*(1-UserAcc[i]))/private$sumfil(self$Values)[i])
         ConfInt[[i]]<-c(private$ConfInt(UserAcc[i],VarUserAcc[i],a)$ConfInt_inf,
                         private$ConfInt(UserAcc[i],VarUserAcc[i])$ConfInt_sup,a)
       }
     return(list(UserAcc=UserAcc,VarUserAcc=VarUserAcc,Conf_Int=ConfInt))
     },




      #' @description Public method for deriving the index called
      #' user’s accuracy for a specific class \eqn{i} in a ConfMatrix object
      #' instance. The user’s accuracy for the class \eqn{i} of a thematic
      #' product is calculated by dividing the value in the diagonal of
      #' class \eqn{i} by the sum of all values in the row of the class i
      #' (row marginal). The method also offers the variance and confidence
      #' interval. The reference \insertCite{congalton2008;textual}{ConfMatrix}
      #' is followed for the computations.
      #' @description
      #'  \deqn{
      #' UserAcc_{i}=\dfrac{x_{ii}}{x_{i+}}
      #' }
      #' \deqn{
      #' \sigma^2_{UserAcc_i}=\dfrac{UserAcc_i \cdot (1-UserAcc_i)}{N_{i}}
      #' }
      #'
      #' @param i \verb{
      #' Class to evaluate, where} \eqn{i \in \mathbb{Z}-\{0\}}.
      #'
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #'}
      #' @return A list of real values containing the user’s accuracy
      #' for class i, its variance, and its confidence interval.
      #'
      #' @examples
      #' A<-matrix(c(65,6,0,4,4,81,11,7,22,5,85,3,24,8,19,90), nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A,Source="Congalton and Green 2008")
      #' p$UserAcc_i(2)
      #'
      #' @aliases NULL

     UserAcc_i=function(i,a=NULL){
      UserAcc_i <- self$Values[i,i] / private$sumfil(self$Values)[i]
      VarUserAcc_i <- abs((UserAcc_i*(1-UserAcc_i))/private$sumfil(self$Values)[i])
      ConfInt <- private$ConfInt(UserAcc_i,VarUserAcc_i,a)
     return(list(UserAcc_i=UserAcc_i,VarUserAcc_i=VarUserAcc_i,
                 Conf_Int=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup)))
     },



      #' @description Public method that provides the arithmetic average,
      #' without weighing, of all user’s accuracies of a ConfMatrix object
      #' instance. The method also offers the variance and confidence
      #' interval. The reference \insertCite{tung1988;textual}{ConfMatrix} is
      #' followed for the calculations.
      #' @description
      #'  \deqn{
      #' AvUserAcc=\dfrac{1}{M} \sum^M_{i=1} \dfrac{x_{ii}}
      #' { x_{i+}}
      #' }
      #'\deqn{
      #' \sigma^2_{AvUserAcc}=\dfrac{AvUserAcc \cdot (1-AvUserAcc)}{N_{Total}}
      #' }
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #' @return
      #' A list of real values containing the average
      #' user’s accuracy, its variance, and its confidence interval.
      #'
      #' @examples
      #' A<-matrix(c(352,43,89,203),nrow=2,ncol=2)
      #' p<-ConfMatrix$new(A,Source="Tung and LeDrew 1988")
      #' p$AvUserAcc()
      #'
      #' @aliases NULL

     AvUserAcc = function(a=NULL){
       for (i in 1:length(private$sumfil(self$Values))) {
          if (private$sumfil(self$Values)[i] == 0) {
            stop ("/ by 0")
          }
       }
      AvUserAcc <- 1/sqrt(length(self$Values)) *
        sum (diag(self$Values)/private$sumfil(self$Values))
      VarAvUserAcc <- abs((AvUserAcc*(1-AvUserAcc))/sum(self$Values))
      ConfInt <- private$ConfInt(AvUserAcc,VarAvUserAcc,a)
     return(list(AvUserAcc=AvUserAcc,VarAvUserAcc=VarAvUserAcc,
            Conf_Int=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup)))
     },


      #' @description Public method that provides the combined user's accuracy.
      #' Which is the average of the overall accuracy and the average user's
      #' accuracy. The method also offers the
      #' variance and confidence interval. The reference
      #' \insertCite{tung1988;textual}{ConfMatrix} is followed for the calculations.
      #' @description
      #'  \deqn{
      #' CombUserAcc=\dfrac{OverallAcc+AvUserAcc}{2}
      #' }
      #' \deqn{
      #' \sigma^2_{CombUserAcc}=\dfrac{CombUserAcc \cdot (1-CombUserAcc)}{N_{Total}}
      #' }
      #'
      #' where:
      #' \enumerate{
      #'   \item \eqn{OverallAcc}: overall accuracy.
      #'   }
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #' @return A list of real values containing the combined
      #' accuracy from the user's perspective,
      #' its variation and confidence interval.
      #'
      #' @examples
      #' A<-matrix(c(352,43,89,203),nrow=2,ncol=2)
      #' p<-ConfMatrix$new(A,Source="Tung and LeDrew 1988")
      #' p$CombUserAcc()
      #'
      #' @aliases NULL

     CombUserAcc = function(a=NULL){
      CombUserAcc <- (self$OverallAcc()[[1]] + self$AvUserAcc()[[1]]) / 2
      VarCombUserAcc <- abs((CombUserAcc*(1-CombUserAcc))/sum(self$Values))
      ConfInt <- private$ConfInt(CombUserAcc,VarCombUserAcc,a)
     return(list(CombUserAcc=CombUserAcc,VarCombUserAcc=VarCombUserAcc,
            Conf_Int=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup)))
     },


# ProdAcc -----------------------------------------------------------------


      #' @description Public method for deriving the index called
      #' producer’s accuracy for all the classes in a ConfMatrix
      #' object instance. The producer’s accuracy for the class i
      #' of a thematic product is calculated by dividing the value
      #' in the diagonal of class \eqn{i} by the sum of all values in the
      #' row of the class \eqn{i} (column marginal). The method also
      #' offers the variance and confidence interval. The reference
      #' \insertCite{congalton2008;textual}{ConfMatrix} if followed for the
      #' computations.
      #' @description
      #'  \deqn{
      #' ProdAcc=\dfrac{x_{ii}}{x_{+j}}
      #' }
      #'\deqn{
      #' \sigma^2_{ProdAcc}=\dfrac{ProdAcc \cdot (1-ProdAcc)}{N_{j}}
      #' }
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #' @return A list of vectors each one containing the producer’s
      #' accuracy real values for all classes, their variances and
      #' confidence intervals for each class, respectively.
      #'
      #' @examples
      #' A<-matrix(c(65,6,0,4,4,81,11,7,22,5,85,3,24,8,19,90),nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A,Source="Congalton and Green 2008")
      #' p$ProdAcc()
      #'
      #' @aliases NULL

     ProdAcc = function (a=NULL){
      n <- sqrt(length(self$Values))
      ProdAcc <- rep(0,n)
      VarProdAcc<-rep(0,n)
      ConfInt<-list()
        for(i in 1:n){
          ProdAcc[i] <- self$Values[i,i] / private$sumcol(self$Values)[i]
          VarProdAcc[i]<-abs((ProdAcc[i]*(1-ProdAcc[i]))/private$sumcol(self$Values)[i])
          ConfInt[[i]]<-c(private$ConfInt(ProdAcc[i],VarProdAcc[i],a)$ConfInt_inf,
                          private$ConfInt(ProdAcc[i],VarProdAcc[i])$ConfInt_sup,a)
        }
     return(list(ProdAcc=ProdAcc,VarProdAcc=VarProdAcc,Conf_Int=ConfInt))
     },




      #' @description Public method for deriving the index called
      #' producer’s accuracy for a specific class \eqn{i} in a ConfMatrix
      #' object instance. The user’s accuracy for the class \eqn{i} of a
      #' thematic product is calculated by dividing the value in the
      #' diagonal of class \eqn{i} by the sum of all values in the column
      #' of the class \eqn{i} (column marginal). The method also offers
      #' the variance and confidence interval. The reference
      #' \insertCite{congalton2008;textual}{ConfMatrix} is followed for the
      #' calculations.
      #' @description
      #'  \deqn{
      #' ProdAcc_{i}=\dfrac{x_{ii}}{x_{+j}}
      #' }
      #'\deqn{
      #' \sigma^2_{ProdAcc_i}=\dfrac{ProdAcc_i \cdot (1-ProdAcc_i)}{N_{j}}
      #' }
      #'
      #' @param i \verb{
      #' Producer class to evaluate, where} \eqn{i \in \mathbb{Z}-\{0\}.
      #' }
      #'
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #'
      #' @return A list of real values containing the producer’s
      #' accuracy for class i, its variance, and its confidence
      #' interval.
      #'
      #' @examples
      #' A<-matrix(c(65,6,0,4,4,81,11,7,22,5,85,3,24,8,19,90),nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A,Source="Congalton and Green 2008")
      #' p$ProdAcc_i(1)
      #'
      #' @aliases NULL

     ProdAcc_i = function(i,a=NULL){
      ProdAcc_i <- self$Values[i,i] / private$sumcol(self$Values)[i]
      VarProdAcc_i <- abs((ProdAcc_i*(1-ProdAcc_i))/private$sumcol(self$Values)[i])
      ConfInt <- private$ConfInt(ProdAcc_i,VarProdAcc_i,a)
     return(list(ProdAcc_i=ProdAcc_i,VarProdAcc_i=VarProdAcc_i,
                 Conf_Int=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup)))
     },




      #' @description Public method that provides the arithmetic
      #' average of all producer’s accuracies of a ConfMatrix object
      #' instance. The method also offers the variance and confidence
      #' interval. The reference \insertCite{tung1988;textual}{ConfMatrix}
      #' is followed for the calculations.
      #' @description
      #'  \deqn{
      #' AvProdAcc=\dfrac{1}{M} \sum^M_{i=1} \dfrac{x_{ii}}
      #' { x_{+j}}
      #' }
      #'\deqn{
      #' \sigma^2_{AvProdAcc}=\dfrac{AvProdAcc \cdot (1-AvProdAcc)}{N_{Total}}
      #' }
      #'
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #' @return A list of real values containing the average producer’s
      #' accuracy, its variance, and its confidence interval.
      #'
      #' @examples
      #' A<-matrix(c(352,43,89,203),nrow=2,ncol=2)
      #' p<-ConfMatrix$new(A,Source="Tung and LeDrew 1988")
      #' p$AvProdAcc()
      #'
      #' @aliases NULL

     AvProdAcc = function(a=NULL){
      AvProdAcc <- 1/sqrt(length(self$Values)) *
        sum (diag(self$Values)/private$sumcol(self$Values))
      VarAvProdAcc <- abs((AvProdAcc*(1-AvProdAcc))/sum(self$Values))
      ConfInt <- private$ConfInt(AvProdAcc,VarAvProdAcc,a)
     return(list(AvProdAcc=AvProdAcc,VarAvProdAcc=VarAvProdAcc,
                 Conf_Int=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup)))
     },




      #' @description Public method that provides the combined producer's
      #' accuracy. Which is the average of the overall accuracy and the average
      #' producer accuracy. The method also offers the
      #' variance and confidence interval. The reference
      #' \insertCite{tung1988;textual}{ConfMatrix} is followed for the calculations.
      #' @description
      #'  \deqn{
      #' CombProdAcc=\dfrac{OverallAcc+AvProdAcc}{2}
      #' }
      #'\deqn{
      #' \sigma^2_{CombProdAcc}=\dfrac{CombProdAcc \cdot (1-CombProdAcc)}{N_{Total}}
      #' }
      #' where:
      #'
      #' \enumerate{
      #'   \item \eqn{OverallAcc}: overall accuracy.
      #'   \item \eqn{AvProdAcc}: average accuracy from producer's perspective.
      #' }
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #' @return A list of real values containing the combined accuracy
      #' from producer's perspective, its variance and confidence interval.
      #'
      #' @examples
      #' A<-matrix(c(352,43,89,203),nrow=2,ncol=2)
      #' p<-ConfMatrix$new(A,Source="Tung and LeDrew 1988")
      #' p$CombProdAcc()
      #'
      #' @aliases NULL

     CombProdAcc = function(a=NULL){
      CombProdAcc <- (self$OverallAcc()[[1]] + self$AvProdAcc()[[1]]) / 2
      VarCombProdAcc <- abs((CombProdAcc*(1-CombProdAcc))/sum(self$Values))
      ConfInt <- private$ConfInt(CombProdAcc,VarCombProdAcc,a)
     return(list(CombProdAcc=CombProdAcc,VarCombProdAcc=VarCombProdAcc,
            Conf_Int=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup)))
     },


# UserProd ----------------------------------------------------------------




      #' @description Public method that calculates the user’s and the
      #' producer’s indexes jointly. This method is equivalent to the methods
      #' ConfMatrix$UserAcc and ConfMatrix$ProdAcc.
      #' @return A list containing the producer's and user's accuracies and
      #' their standard deviations, respectively.
      #' @examples
      #' A<-matrix(c(65,6,0,4,4,81,11,7,22,5,85,3,24,8,19,90),nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A,Source="Congalton and Green 2008")
      #' p$UserProdAcc()
      #'
      #' @aliases NULL

     UserProdAcc =function(){
      nc <- nrow(self$Values)
       for (i in 1:nc){
         pcpa <- self$ProdAcc()[[1]]
         pcua <-self$UserAcc()[[1]]
         pcpasd <- sqrt(self$ProdAcc()[[2]])
         pcuasd <- sqrt(self$UserAcc()[[2]])
       }
     return(list(ProdAcc=pcpa,UserAcc= pcua,ProdAccSDeviation=pcpasd,
                 UserAccSDeviation=pcuasd))
     },



      #' @description Public method that provides the combined accuracy,
      #' defined by the average of the overall accuracy and the Hellden's
      #' average accuracy, which refers to the average user's and producer's
      #' accuracies. The method also offers the
      #' variance and confidence interval. The reference
      #' \insertCite{liu2007;textual}{ConfMatrix} is followed for the calculations.
      #' @description
      #'  \deqn{
      #' CombUserProdAcc=\dfrac{OverallAcc+AvHellAcc}{2}
      #' }
      #'  \deqn{
      #' \sigma^2_{CombUserProdAcc}=\dfrac{CombUserProdAcc \cdot
      #' (1-CombUserProdAcc)}{N_{Total}}
      #' }
      #' where:
      #'
      #' \enumerate{
      #'   \item \eqn{OverallAcc}: overall accuracy.
      #'   \item \eqn{AvHellAcc}: average of Hellden's mean accuracy index.
      #' }
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #' @return A list of real values containing the combined accuracy from both user's
      #' and producer's perspectives, its variance and confidence interval.
      #'
      #' @examples
      #' A<-matrix(c(65,6,0,4,4,81,11,7,22,5,85,3,24,8,19,90),nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A,Source="Congalton and Green 2008")
      #' p$CombUserProdAcc()
      #'
      #' @aliases NULL

     CombUserProdAcc = function(a=NULL){
      CombUserProdAcc <- (self$OverallAcc()[[1]]+self$AvHellAcc()[[1]])/2
      VarCombUserProdAcc <- abs((CombUserProdAcc*
                            (1-CombUserProdAcc))/sum(self$Values))
      ConfInt <- private$ConfInt(CombUserProdAcc,VarCombUserProdAcc,a)
     return(list(CombUserProdAcc=CombUserProdAcc,
                 VarCombUserProdAcc=VarCombUserProdAcc,
                 Conf_Int=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup)))
     },



      #' @description Public method that provides the arithmetic
      #' average of all user’s and producer’s accuracy indexes of
      #' a ConfMatrix object instance. The method also offers the
      #' variance and confidence interval. The reference
      #' \insertCite{liu2007;textual}{ConfMatrix} is followed for the
      #' calculations.
      #' @description
      #'  \deqn{
      #' AvUserProdAcc=\dfrac{AvUserAcc+AvProdAcc}{2}
      #' }
      #'  \deqn{
      #' \sigma^2_{AvUserProdAcc}=\dfrac{AvUserProdAcc \cdot (1-AvUserProdAcc)}{N_{Total}}
      #' }
      #'
      #' where:
      #'
      #' \enumerate{
      #'   \item \eqn{AvUserAcc}: average accuracy from user's perspective.
      #'   \item \eqn{AvProdAcc}: average accuracy from producer's perspective.
      #' }
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #' @return A list of real values containing the average mean
      #' precision values from the user's and producer's perspective,
      #' their variance and confidence interval.
      #'
      #' @examples
      #' A<-matrix(c(65,6,0,4,4,81,11,7,22,5,85,3,24,8,19,90),nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A,Source="Congalton and Green 2008")
      #' p$AvUserProdAcc()
      #'
      #' @aliases NULL

     AvUserProdAcc = function(a=NULL){
      AvUserProdAcc <- (self$AvUserAcc()[[1]] + self$AvProdAcc()[[1]]) / 2
      VarAvUserProdAcc <- abs((AvUserProdAcc*(1-AvUserProdAcc))/sum(self$Values))
      ConfInt <- private$ConfInt(AvUserProdAcc,VarAvUserProdAcc,a)
     return(list(AvUserProdAcc=AvUserProdAcc,VarAvUserProdAcc=VarAvUserProdAcc,
                 Conf_Int=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup)))
     },



      #' @description Public method that provides the average of
      #' user’s and producer’s accuracies for a specific class i.
      #' The method also offers the variance and confidence
      #' interval. The reference \insertCite{liu2007;textual}{ConfMatrix}
      #' is followed for the calculations.
      #'  \deqn{
      #' AvUserProdAcc_i=\dfrac{UserAcc_i+ProdAcc_i}{2}
      #' }
      #'\deqn{
      #' \sigma^2_{AvUserProdAcc_i}=\dfrac{AvUserProdAcc_i
      #' \cdot (1-AvUserProdAcc_i)}{N_{ij}}
      #' }
      #' where:
      #'
      #' \enumerate{
      #'   \item \eqn{UserAcc_i}: user accuracy index for class i.
      #'   \item \eqn{ProdAcc_i}: producer accuracy index for class i.
      #' }
      #'
      #' @param i \verb{
      #' Class to evaluate, where} \eqn{i \in \mathbb{Z}-\{0\}}.
      #'
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #'
      #' @return A list of real values containing the average of
      #' user’s and producer’s accuracies, its variance and
      #' confidence interval for class i.
      #'
      #' @examples
      #' A<-matrix(c(65,6,0,4,4,81,11,7,22,5,85,3,24,8,19,90),nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A,Source="Congalton and Green 2008")
      #' p$AvUserProdAcc_i(2)
      #'
      #' @aliases NULL

     AvUserProdAcc_i = function(i,a=NULL){
      AvUserProdAcc_i <- (self$UserAcc_i(i)[[1]] + self$ProdAcc_i(i)[[1]])/2
      VarAvUserProdAcc_i <- abs((AvUserProdAcc_i*
                          (1-AvUserProdAcc_i))/(private$sumcol(self$Values)[i]+private$sumfil(self$Values)[i]-self$Values[i,i]))
      ConfInt <- private$ConfInt(AvUserProdAcc_i,VarAvUserProdAcc_i,a)
     return(list(AvUserProdAcc_i=AvUserProdAcc_i,
                 VarAvUserProdAcc_i=VarAvUserProdAcc_i,
                 ConfInt=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup)))
     },



      #' @description Public method that calculates the weighted user's,
      #' producer’s and overall accuracies and their standard deviations.
      #' The reference \insertCite{congalton2008;textual}{ConfMatrix} is followed
      #' for the computations.
      #'
      #' Be
      #' \deqn{
      #' Overall_W=\dfrac{\sum^M_{i=1} p_{ii} }{\sum^M_{i,j=1} p_{ij}}
      #' }
      #
      #' where \eqn{p_{ij}=\dfrac{x_{ij}}{\sum^M_{i,j=1} x_{ij}}}
      #
      #' \deqn{
      #' UserAcc_W=\dfrac{p_{o_{i+}}}{p_{i+}}
      #' }
      #' \deqn{
      #' ProdAcc_W=\dfrac{p_{o_{+j}}}{p_{+j}}
      #' }
      #'
      #' \deqn{
      #' \sigma^2_{UserAcc_W}=\dfrac{UserAcc_W \cdot (1-UserAcc_W)}{N_i}
      #' }
      #' \deqn{
      #' \sigma^2_{ProdAcc_W}=\dfrac{ProdAcc_W \cdot (1-ProdAcc_W)}{N_j}
      #' }
      #'
      #'
      #' where \eqn{p_o=\sum^M_{i,j=1} w_{ij}p_{ij}} and \eqn{0 \leq w_{ij} \leq 1}
      #' for \eqn{i \neq j} and \eqn{w_{ii}=1} for \eqn{i = j}.
      #'
      #'
      #' @param WM Weight matrix (as matrix)
      #' @return A list with the weight matrix, the product of the
      #' confusion matrix and the weight matrix, overall, user and
      #' producer weighted accuracies and their standard deviations.
      #' @examples
      #' A<-matrix(c(65,6,0,4,4,81,11,7,22,5,85,3,24,8,19,90),nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A,Source="Congalton and Green 2008")
      #' WM<- t(matrix(c(1,0,0.67,1,0,1,0,0,1,0,1,1,0.91,0,0.61,1),nrow=4,ncol=4))
      #' p$UserProdAcc_W(WM)
      #'
      #' @aliases NULL

     UserProdAcc_W =function(WM){
        ncol <- private$sumcol(self$Values)
        nrow<- private$sumfil(self$Values)
        ConfM<- self$Values/sum(self$Values)

        if(any(WM)>1){
          WM<-WM/sum(WM)
        }
        diag(WM)<-1
        WConfM<-ConfM*WM
        WOverallAcc <- sum(diag(WConfM))/sum(WConfM)
        mcol<- apply(WConfM,2,sum)
        mrow<- apply(WConfM,1,sum)
        pj <- apply(ConfM,2,sum)
        pi <- apply(ConfM,1,sum)

        nc <- nrow(ConfM)
        wi<- matrix(rep(0, nc), nrow =1, ncol=nc, byrow=TRUE)
        wj<- matrix(rep(0, nc), nrow =1, ncol=nc, byrow=TRUE)
        wpcua<- matrix(rep(0, nc), nrow =1, ncol=nc, byrow=TRUE)
        wpcpa<- matrix(rep(0, nc), nrow =1, ncol=nc, byrow=TRUE)
        pcpasd <-matrix(rep(0, nc), nrow =1, ncol=nc, byrow=TRUE)
        pcuasd <-matrix(rep(0, nc), nrow =1, ncol=nc, byrow=TRUE)

            for (i in 1:nc){
                wi[i]    <- sum(pj*WM[i,])
                wj[i]    <- sum(pi*WM[,i])
                wpcua[i]  <- sum(WConfM[i,])/sum(ConfM[i,])
                wpcpa[i]  <- sum(WConfM[,i])/sum(ConfM[,i])
                pcpasd[i] <- sqrt(wpcpa[i]*(1-wpcpa[i])/ncol[i])
                pcuasd[i] <- sqrt(wpcua[i]*(1-wpcua[i])/nrow[i])
            }

     return(list(OriginalWeightMatrix=WM,WMatrix=WConfM, WOverallAcc=WOverallAcc,
                 WPrAcc=wpcpa,WPrAccSDeviation=pcpasd,WUserAcc= wpcua,
                 WUserAccSDeviation=pcuasd))
     },


# Functions that return indices, variances and confidence interval --------


      #' @description Public method to calculate the global index called
      #' Overall Accuracy. The Overall Accuracy is calculated by dividing
      #' the sum of the entries that form the major diagonal (i.e., the
      #' number of correct classifications) by the total number of cases.
      #' The method also offers the variance and confidence interval.
      #' The reference \insertCite{congalton2008;textual}{ConfMatrix}
      #' is followed for the computations.
      #'
      #' \deqn{
      #' OverallAcc = \dfrac{\sum_{i=1}^{M} x_{ii}}{\sum_{i, j=1}^{M} x_{ij}}
      #' }
      #'
      #' \deqn{
      #' \sigma^2_{OverallAcc}=\dfrac{OverallAcc \cdot (1-OverallAcc)}{N_{Total}}
      #' }
      #'
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #' @return A list of real values containing the overall accuracy,
      #' its variance, and its confidence interval.
      #'
      #'
      #'
      #' @examples
      #' A<-matrix(c(65,6,0,4,4,81,11,7,22,5,85,3,24,8,19,90),nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A)
      #' p$OverallAcc()
      #'
      #' @aliases NULL

     OverallAcc = function(a=NULL) {
     index <- sum(diag(self$Values))/sum(self$Values)
     VarIndex<-abs((index*(1-index))/sum(self$Values))
     ConfInt<-private$ConfInt(index,VarIndex,a)
     return(list(OverallAcc=index,VarOverallAcc=VarIndex,
                 Conf_Int=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup)))
     },


# Kappa -------------------------------------------------------------------



      #' @description Public method that provides Kappa coefficient,
      #' which measures the relationship between the observed proportion
      #' of agreement and the proportion expected to occur by chance.
      #' The method also offers the variance and confidence interval.
      #' The reference \insertCite{cohen1960;textual}{ConfMatrix} is followed
      #' for the calculations.
      #' @description
      #' \deqn{
      #' Kappa=\dfrac{OverallAcc-ExpAcc}{1-ExpAcc}
      #' }
      #' \deqn{
      #' ExpAcc= \dfrac{x_{+ i}x_{i +}}{ ( \sum_{i,j=1}^M x_{ij} )^{2} }
      #' }
      #' \deqn{
      #' \sigma^2_{Kappa}=\dfrac{OverallAcc-ExpAcc}{(1-ExpAcc)^2 \cdot N_{Total}}
      #' }
      #'
      #'
      #'
      #' where:
      #'
      #' \enumerate{
      #'   \item \eqn{OverallAcc}: overall accuracy.
      #'   \item \eqn{ExpAcc}: expected accuracy of agreement if agreement
      #'   were purely random.
      #' }
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #' @return A list of real values containing with kappa
      #' coefficient, its variance and confidence interval.
      #'
      #' @examples
      #' A<-matrix(c(65,6,0,4,4,81,11,7,22,5,85,3,24,8,19,90),nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A,Source="Congalton and Green 2008")
      #' p$Kappa()
      #'
      #' @aliases NULL

     Kappa = function(a=NULL){
      ExpAcc <- (sum (private$sumfil(self$Values) * private$sumcol(self$Values)))/sum(self$Values)^2
      if (1-ExpAcc == 0){
       stop ("/ by 0")
      }else{
        kappa <- (self$OverallAcc()[[1]]- ExpAcc) / (1 - ExpAcc)
        VarKappa <- abs((self$OverallAcc()[[1]]*(1-self$OverallAcc()[[1]]))
                    /(sum(self$Values)*(1-ExpAcc)^2))
        ConfInt <- private$ConfInt(kappa,VarKappa,a)
      }
     return(list(Kappa=kappa,VarKappa=VarKappa,
            Conf_Int=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup)))
     },


      #' @description Public method that provides the overall modified
      #' kappa coefficient. The method also offers the
      #' variance and confidence interval. The references
      #' \insertCite{stehman1997;textual}{ConfMatrix} and \insertCite{foody1992;textual}{ConfMatrix}
      #' are followed for the calculations.
      #'
      #'  \deqn{
      #' ModKappa=\dfrac{OverallAcc-\dfrac{1}{M}}{1-\dfrac{1}{M}}
      #' }
      #' \deqn{
      #' \sigma^2_{ModKappa}=\dfrac{OverallAcc \cdot (1- OverallAcc)}
      #' { \left(1-\dfrac{1}{M} \right)^2 \cdot N_{Total}}
      #' }
      #' where:
      #'
      #' \enumerate{
      #'   \item \eqn{OverallAcc}: overall accuracy.
      #' }
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #' @return A list of real values containing modified coefficient
      #' kappa, its variance and its confidence interval.
      #' @examples
      #' A <- matrix(c(317,61,2,35,23,120,4,29,0,0,60,0,0,0,0,8),nrow=4,ncol=4)
      #' p <- ConfMatrix$new(A,Source="Foody 1992")
      #' p$ModKappa()
      #'
      #' @aliases NULL

     ModKappa = function(a=NULL){
       ModKappa <- (self$OverallAcc()[[1]] -
                   1/sqrt(length(self$Values)))/(1 -
                   1/sqrt(length(self$Values)))
       VarModKappa <- (self$OverallAcc()[[1]]*(1-self$OverallAcc()[[1]]))/
         (((1 - 1/sqrt(length(self$Values)))^2)*sum(self$Values))
       ConfInt <- private$ConfInt(ModKappa,VarModKappa,a)
     return(list(ModKappa=ModKappa,VarModKappa=VarModKappa,
                 Conf_Int=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup)))
     },



      #' @description Public method derived by the kappa coefficient evaluated
      #' from the user's perspective, for a specific class i. The method
      #' also offers the variance and confidence interval. The reference
      #' \insertCite{rosenfield1986;textual}{ConfMatrix} is followed
      #' for the calculations.
      #' @description
      #'  \deqn{
      #' UserKappa_i=\dfrac{UserAcc_i-\dfrac{ x_{i + }}
      #' {\sum^M_{i,j=1} x_{ij}}}{1-\dfrac{ x_{i + }}
      #' {\sum^M_{i,j=1} x_{ij}}}
      #' }
      #'
      #'  \deqn{
      #' \sigma^2_{UserKappa_i}=\dfrac{UserAcc_i \cdot (1-UserAcc_i)}
      #' { \left(1-\dfrac{ x_{i + }}
      #' {\sum^M_{i,j=1} x_{ij}}\right)^2 \cdot N_{i}}
      #' }
      #' where:
      #'
      #' \enumerate{
      #'   \item \eqn{UserAcc_i}: user accuracy index for class i.
      #' }
      #'
      #' @param i \verb{
      #' Class to evaluate, where} \eqn{i \in \mathbb{Z}-\{0\}}.
      #'
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #' @return A list of real values containing the kappa coefficient
      #' (user’s perspective), its variance and its confidence interval.
      #'
      #' @examples
      #' A<-matrix(c(73,13,5,1,0,21,32,13,3,0,16,39,35, 29,13,3,5,7,28,48,1,0,2,3,17),
      #' nrow=5,ncol=5)
      #' p<-ConfMatrix$new(A,Source="Næsset 1996")
      #' p$UserKappa_i(2)
      #'
      #' @aliases NULL

     UserKappa_i = function(i,a=NULL){
      if (1 - private$sumcol(self$Values)[i]/sum(self$Values) == 0) {
       stop ("/ by 0")
      }else{
        UserKappa_i <- (self$UserAcc_i(i)[[1]] -
                        private$sumcol(self$Values)[i]/sum(self$Values)) /
                        (1 - private$sumcol(self$Values)[i]/sum(self$Values))
        VarUserKappa_i <- (self$UserAcc_i(i)[[1]]*(1-self$UserAcc_i(i)[[1]])) /
          (((1 - private$sumcol(self$Values)[i]/sum(self$Values))^2)*private$sumfil(self$Values)[i])
        ConfInt <- private$ConfInt(UserKappa_i,VarUserKappa_i,a)
        }

     return(list(UserKappa_i=UserKappa_i,VarUserKappa_i=VarUserKappa_i,
                 Conf_Int=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup)))
     },



      #' @description Public method, derived from the general modified
      #' kappa coefficient, which provides the modified kappa coefficient
      #' from the user's perspective and for a specific class i. Equitable probabilities
      #' of belonging to each class are assumed. The method also offers
      #' the variance and confidence interval. The references
      #' \insertCite{stehman1997;textual}{ConfMatrix} and \insertCite{foody1992;textual}{ConfMatrix}
      #' are followed for the calculations.
      #'  \deqn{
      #' ModKappaUser_i=\dfrac{UserAcc_i-\dfrac{1}{M}}
      #' {1-\dfrac{1}{M}}
      #' }
      #' \deqn{
      #' \sigma^2_{ModKappaUser_i}=\dfrac{UserAcc_i
      #' \cdot (1- UserAcc_i)}{ \left(1- \dfrac{1}{M} \right)^2 \cdot N_{i}}
      #' }
      #' where:
      #'
      #' \enumerate{
      #'   \item \eqn{UserAcc_i}: user accuracy index for class i.
      #' }
      #' @param i \verb{
      #' Class to evaluate, where} \eqn{i \in \mathbb{Z}-\{0\}}.
      #'
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #' @return A list of real values containing the modified kappa
      #' coefficient from the user's perspective, its variance and
      #' confidence interval.
      #'
      #' @examples
      #' A<-matrix(c(0,12,0,0,12,0,0,0,0,0,0,12,0,0,12,0),nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A,Source="Liu et al. 2007")
      #' p$ModKappaUser_i(2)
      #'
      #' @aliases NULL

     ModKappaUser_i = function(i,a=NULL){
       ModKappaUser_i <- (self$UserAcc_i(i)[[1]] -
                         1/sqrt(length(self$Values)))/(1 -
                         1/sqrt(length(self$Values)))
       VarModKappaUser_i <- (self$UserAcc_i(i)[[1]]*(1-self$UserAcc_i(i)[[1]]))/
         (((1 - 1/sqrt(length(self$Values)))^2)*(private$sumfil(self$Values)[i]))
       ConfInt <- private$ConfInt(ModKappaUser_i,VarModKappaUser_i,a)
     return(list(ModKappaUser_i=ModKappaUser_i,
                 VarModKappaUser_i=VarModKappaUser_i,
                 Conf_Int=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup)))
     },



      #' @description Public method derived by the kappa coefficient evaluated
      #' from the producer's perspective, for a specific class i. The method
      #' also offers the variance and confidence interval. The reference
      #' \insertCite{rosenfield1986;textual}{ConfMatrix} is followed
      #' for the calculations.
      #' @description
      #'  \deqn{
      #' ProdKappa_i=\dfrac{ProdAcc_i-\dfrac{ x_{ + i }}
      #' {\sum^M_{i,j=1} x_{ij}}}{1-\dfrac{ x_{+ i }}
      #' {\sum^M_{i,j=1} x_{ij}}}
      #' }
      #'
      #'  \deqn{
      #' \sigma^2_{ProdKappa_i}=\dfrac{ProdAcc_i \cdot (1- ProdAcc_i)}
      #' {\left(1-\dfrac{ x_{+ i }}
      #' {\sum^M_{i,j=1} x_{ij}} \right)^2 \cdot N_{j}}
      #' }
      #' where:
      #'
      #' \enumerate{
      #'   \item \eqn{ProdAcc_i}: producer accuracy index for class i.
      #' }
      #'
      #' @param i \verb{
      #' Class to evaluate, where} \eqn{i \in \mathbb{Z}-\{0\}}.
      #'
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #' @return A list of real values containing the coefficient kappa
      #' (producer’s), its variance and its confidence interval.
      #' @examples
      #' A <- matrix(c(73,13,5,1,0,21,32,13,3,0,16,39,35,29,13,3,5,7,28,48,1,0,2,3,17),
      #' nrow=5,ncol=5)
      #' p<-ConfMatrix$new(A,Source="Næsset 1996")
      #' p$ProdKappa_i(2)
      #'
      #' @aliases NULL

      ProdKappa_i = function(i,a=NULL){
      if (1 - private$sumfil(self$Values)[i]/sum(self$Values) == 0) {
       stop ("/ by 0")
      }else{
        ProdKappa_i <- (self$ProdAcc_i(i)[[1]] - private$sumfil(self$Values)[i]/sum(self$Values)) / (1 - private$sumfil(self$Values)[i]/sum(self$Values))
        VarProdKappa_i <-(self$ProdAcc_i(i)[[1]]*(1-self$ProdAcc_i(i)[[1]])) /
          (((1 - private$sumfil(self$Values)[i]/sum(self$Values))^2)*private$sumcol(self$Values)[i])
        ConfInt <- private$ConfInt(ProdKappa_i,VarProdKappa_i,a)
        }
     return(list(ProdKappa_i=ProdKappa_i,VarProdKappa_i=VarProdKappa_i,
                 Conf_Int=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup)))
     },



      #' @description Public method, derived from the general modified
      #' kappa coefficient, which provides the modified kappa coefficient
      #' from the producer's perspective and for a specific class i. Equitable
      #' probabilities of belonging to each class are assumed. The method also
      #' offers the variance and confidence interval. The references
      #' \insertCite{stehman1997;textual}{ConfMatrix} and \insertCite{foody1992;textual}{ConfMatrix}
      #' are followed for the calculations.
      #' @description
      #'  \deqn{
      #' ModKappaProd_i=\dfrac{ProdAcc_i-\dfrac{1}{M}}
      #' {1-\dfrac{1}{M}}
      #' }
      #' \deqn{
      #' \sigma^2_{ModKappaProd_i}=\dfrac{ProdAcc_i
      #' \cdot (1- ProdAcc_i)}{ \left( 1-\dfrac{1}{M} \right)^2 \cdot N_{j}}
      #' }
      #' where:
      #'
      #' \enumerate{
      #'   \item \eqn{ProdAcc_i}: producer accuracy index for class i.
      #' }
      #' @param i \verb{
      #' Class to evaluate, where} \eqn{i \in \mathbb{Z}-\{0\}}.
      #'
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #' @return
      #' A list of real values containing the modified kappa coefficient
      #' from the producer's perspective, its variance and confidence interval.
      #'
      #' @examples
      #' A<-matrix(c(317,61,2,35,23,120,4,29,0,0,60,0,0,0,0,8),nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A,Source="Foody 1992")
      #' p$ModKappaProd_i(2)
      #'
      #' @aliases NULL

     ModKappaProd_i = function(i,a=NULL){
      ModKappaProd_i <- (self$ProdAcc_i(i)[[1]] - 1/sqrt(length(self$Values))) / (1 - 1/sqrt(length(self$Values)))
      VarModKappaProd_i <- (self$ProdAcc_i(i)[[1]]*(1-self$ProdAcc_i(i)[[1]]))/
        (((1 - 1/sqrt(length(self$Values)))^2)*(private$sumcol(self$Values)[i]))
      ConfInt <- private$ConfInt(ModKappaProd_i,VarModKappaProd_i,a)
     return(list(ModKappaProd_i=ModKappaProd_i,
                 VarModKappaProd_i=VarModKappaProd_i,
                 Conf_Int=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup)))
     },



      #' @description  Public method that calculates the general Kappa
      #' agreement index, its standard deviation and the test statistic
      #' to test its significance. The delta method has been used to calculate
      #' the sample variance. The reference
      #' \insertCite{congalton2008;textual}{ConfMatrix} is followed for the computations.
      #'
      #' \deqn{
      #' Kappa=\dfrac{OverallAcc-ExpAcc}{1-ExpAcc}
      #' }
      #'
      #'  \deqn{
      #' ExpAcc= \dfrac{x_{+ j} \cdot x_{i +}}{\sum_{(i,j=1}^M x_{ij})^{2}}
      #' }
      #' \deqn{
      #' \sigma^2_{Kappa} = \dfrac{1}{N_{Total}} \left( \dfrac{\theta_1 (1-\theta_1) }{(1-\theta_2)^2}
      #' + \dfrac{2(1-\theta_1)(2\theta_1\theta_2-\theta_3)}{(1-\theta_2)^3}
      #' + \dfrac{(1-\theta_1)^2(\theta_4-4\theta_2^2)}{(1-\theta_2)^4} \right)
      #' }
      #' where
      #'
      #' \deqn{
      #' \theta_1=OverallAcc= \sum_{i, j=1}^{M} \dfrac{ x_{ii}}{ x_{ij}}
      #' }
      #'
      #' \deqn{
      #' \theta_2=ExpAcc=\sum^M_{i=1}
      #' \left( \dfrac{x_{+ i}}{\sum_{j=1}^M x_{ij}}
      #' \cdot \dfrac{x_{i +}}{\sum_{j=1}^M x_{ij}} \right)
      #' }
      #'
      #'
      #' \deqn{
      #' \theta_3=\sum^M_{i=1} \left( \dfrac{x_{ii} x_{+ i}}{\sum_{j=1}^M x_{ij}}
      #' \cdot \dfrac{x_{ii} x_{i +}}{\sum_{j=1}^M x_{ij}} \right)
      #' }
      #'
      #'
      #' \deqn{
      #' \theta_4=\dfrac{1}{ ( \sum_{i,j=1}^M x_{ij})^3} \sum_{i,j=1}^M x_{ij}
      #' (x_{j+}+x_{+i})^2
      #' }
      #'
      #' \deqn{
      #' Z=\dfrac{Kappa}{\sqrt{\sigma^2_{Kappa}}}
      #' }
      #'
      #'
      #' Where:
      #' \enumerate{
      #'   \item \eqn{ExpAcc}: expected accuracy of agreement if agreement
      #'   were purely random.
      #'   \item \eqn{OverallAcc}: overall accuracy.
      #'   \item \eqn{\theta_1, \theta_2, \theta_3, \theta_4}: real values.
      #'   \item \eqn{Z}: the test statistic.
      #' }
      #'
      #'
      #'
      #' @return A list of real values containing the kappa coefficient,
      #' its standard deviation, and the value of its test statistic.
      #' @examples
      #' A<-matrix(c(65,6,0,4,4,81,11,7,22,5,85,3,24,8,19,90),nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A,Source="Congalton and Green 2008")
      #' p$DetailKappa()
      #'
      #' @aliases NULL

     DetailKappa=function (){
       nc <- nrow(self$Values)
       SumaMatriz <-sum(self$Values)
       O1 <- sum(diag(self$Values))/SumaMatriz 
       O2 <- sum((private$sumcol(self$Values)*private$sumfil(self$Values)))/(SumaMatriz*SumaMatriz) 
       O3 <- sum(diag(self$Values)*(private$sumcol(self$Values)+private$sumfil(self$Values)))/(SumaMatriz*SumaMatriz)
       mintermedia1<- matrix(rep(private$sumcol(self$Values), nc), nrow =nc, ncol=nc, byrow=TRUE)
       mintermedia2<- matrix(rep(private$sumfil(self$Values), nc), nrow =nc, ncol=nc, byrow=FALSE)
       mintermedia3 <-(mintermedia1+mintermedia2)^2
       O4 <- sum(self$Values*(t(mintermedia3)) )/(SumaMatriz*SumaMatriz*SumaMatriz)
       t1 <- (1-O1) 
       t2<- (1-O2) 
       K <- (O1-O2)/t2 
       t3<- 2*O1*O2-O3
       t4<- O4-4*(O2^2)
       t5<- O1*t1/(t2^2)+2*t1*t3/(t2^3)+(t1^2)*t4/(t2^4)
       SdK <- sqrt((1/SumaMatriz)*t5) 
       CV <- K/SdK 

     return(list(K=K, SdK=SdK, CV=CV))
     },



      #' @description Public method that calculates the Kappa class agreement
      #' index (conditional Kappa) from the perspective of user (i) and
      #' producer (j) and its standard deviations. The reference
      #' \insertCite{congalton2008;textual}{ConfMatrix} is followed for the computations.
      #'
      #' \deqn{
      #' CondKappa_{user}=\dfrac{\dfrac{x_{ii}}{x_{i+}}-x_{+j}}{1-x_{+j}}
      #' }
      #'
      #' \deqn{
      #' CondKappa_{producer}=\dfrac{\dfrac{x_{ii}}{x_{+j}}-x_{i+}}{1-x_{i+}}
      #' }
      #'
      #' \deqn{
      #' \sigma^2_{CondKappa_{producer}}=\dfrac{1}{N_{Total}} \cdot
      #' \dfrac{x_{+j}-x_{ii}}{x_{+j}^3 (1-x_{i+})^3} \cdot ((x_{+j}-x_{ii})\cdot
      #' (x_{+j}x_{i+}-x_{ii}) + x_{ii} (1-x_{+j}-x_{i+}+x_{ii})  )
      #' }
      #'
      #' \eqn{\sigma^2_{CondKappa_{user}}} is done in an analogous way by exchanging
      #' \eqn{x_{i+}} to \eqn{x_{+j}}.
      #'
      #'
      #'
      #' @return A list of real values containing conditional Kappa index of the user's and the
      #' producer's, and its corresponding standard deviation.
      #' @examples
      #' A<-matrix(c(0.2361,0.0694,0.1389,0.0556,0.1667,0.0417,0.1111,0,0.1806),
      #' ncol=3,nrow=3)
      #' p<-ConfMatrix$new(A,Source="Czaplewski 1994")
      #' p$DetailCondKappa ()
      #'
      #' @aliases NULL


     DetailCondKappa = function(){
       SumaMatriz <-sum(self$Values)
        ConfM<- self$Values/sum(self$Values)
        pcol <- apply(ConfM,2,sum)
        prow<- apply(ConfM,1,sum)
        nc  <- nrow(ConfM)
        Ki_ <- matrix(rep(0, nc), nrow =1, ncol=nc, byrow=TRUE)
        K_j <- matrix(rep(0, nc), nrow =1, ncol=nc, byrow=TRUE)
        Ki_sd <- matrix(rep(0, nc), nrow =1, ncol=nc, byrow=TRUE)
        K_jsd <- matrix(rep(0, nc), nrow =1, ncol=nc, byrow=TRUE)
          for (i in 1:nc){
            Ki_[i] <- ((ConfM[i,i]/prow[i])-pcol[i])/(1-pcol[i])
            K_j[i] <- ((ConfM[i,i]/pcol[i])-prow[i])/(1-prow[i])
            ti1 <- prow[i]-ConfM[i,i]
            ti2<-  ti1/((prow[i]^3)*(1-pcol[i])^3)
            ti3<-  ti1*(pcol[i]*prow[i]-ConfM[i,i])
            ti4<-  ConfM[i,i]*(1-pcol[i]-prow[i]+ConfM[i,i])
            Ki_sd[i] <- sqrt((1/SumaMatriz)*ti2*(ti3+ti4))
            tj1 <- pcol[i]-ConfM[i,i]
            tj2<-  tj1/((pcol[i]^3)*(1-prow[i])^3)
            tj3<-  tj1*(pcol[i]*prow[i]-ConfM[i,i])
            tj4<-  ConfM[i,i]*(1-pcol[i]-prow[i]+ConfM[i,i])
            K_jsd[i] <- sqrt((1/SumaMatriz)*tj2*(tj3+tj4))
          }
     return(list(UserCondKappa=Ki_, ProdCondKappa=K_j,SD_UserCondKappa=Ki_sd,
                 SD_ProdCondKappa=K_jsd))
     },



      #' @description Public method that calculates the general Kappa agreement
      #' index (weighted) and its standard deviation. The reference
      #' \insertCite{fleiss1969,naesset1996;textual}{ConfMatrix} and \insertCite{congalton2008;textual}{ConfMatrix}
      #' are followed for the computations.
      #'
      #' Be \eqn{p_{ij}=\dfrac{x_{ij}}{\sum^M_{i,j} x_{ij}}} for each element \eqn{i,j} for the matrix
      #' and \eqn{0 \leq w_{ij} \leq 1} for \eqn{i \neq j} and \eqn{w_{ii}=1} for \eqn{i = j}.
      #' If the elements of the weight are greater than 1, their value must be given as a percentage.
      #'
      #' Therefore, let:
      #' \deqn{
      #' p_o=\sum^M_{i,j=1} w_{ij}p_{ij}
      #' }
      #' be the weighted agreement, and
      #'
      #' \deqn{
      #' p_c=\sum^M_{i,j=1} w_{ij}p_{i+}p_{+j}
      #' }
      #'
      #' with \eqn{p_{i+}, p_{+j}} analogous to \eqn{x_{i+}, x_{+j}}.
      #'
      #' Then, the weighted Kappa is defined by
      #'
      #' \deqn{
      #' Kappa_w=\dfrac{p_o-p_c}{1-p_c}
      #' }
      #'
      #' The variance may be estimated by
      #'
      #' \deqn{
      #' \sigma^2_{Kappa_w}=\dfrac{1}{N_{Total} (1-p_c)^4} \left(
      #'  \sum^M_{i,j=1} p_{ij} [ w_{ij} (1-p_c)-(\overline{w}_{i+}+\overline{w}_{+j}) (1-p_o)]^2
      #'  -(p_op_c-2p_c+p_o)^2 \right)
      #' }
      #'
      #' where \eqn{\overline{w}_{i+}=\sum^M_{j=1} w_{ij}p_{+j}} and
      #' \eqn{\overline{w}_{+j}=\sum^M_{i=1} w_{ij}p_{i+}}
      #'
      #' Its statistic is given by:
      #'
      #' \deqn{
      #' Z=\dfrac{Kappa_W}{\sqrt{\sigma^2_{Kappa_W}}}
      #' }
      #' @param WM  Weight matrix (as matrix).
      #' @return A list with the weight matrix, kappa index obtained from
      #' the original matrix and the weight matrix, its standard deviations
      #' and the value of its test statistic.
      #' @examples
      #' A <- matrix(c(1,1,0,0,0,5,55,27,23,0,3,30,68,74,4,0,8,8,39,26,0,0,2,4,26),
      #' nrow=5)
      #' WM <- matrix(c(1,0.75,0.5,0.25,0,0.75,1,0.75,0.5,0.25,0.5,0.75,1,0.75,0.5,
      #' 0.25,0.5,0.75,1,0.75,0,0.25,0.5,0.75,1),nrow=5)
      #' p<-ConfMatrix$new(A, Source="Næsset 1996")
      #' p$DetailWKappa(WM)
      #'
      #' @aliases NULL

     DetailWKappa = function(WM){
       nc <- nrow(self$Values)
       SumaMatriz <-sum(self$Values)
       ConfM<-self$Values

       if(any(ConfM)>1){
         ConfM<- self$Values/SumaMatriz
       }
       pcol <- apply(ConfM,2,sum)
       prow<- apply(ConfM,1,sum)
       
       if(any(WM)>1){
        WM<-WM/sum(WM)
        }
        diag(WM)<-1

       Ow1 <- sum(WM*ConfM)
       Ow2 <- sum(t(WM*prow)*pcol)
       c1<- (1-Ow1)
       c2<- (1-Ow2)
       wi_ <- WM %*% pcol
       w_j <- WM %*% prow
       mintermedia1<- matrix(rep(wi_, nc), nrow =nc, ncol=nc, byrow=FALSE)
       mintermedia2<- matrix(rep(w_j, nc), nrow =nc, ncol=nc, byrow=TRUE)
       mintermedia3 <-(mintermedia1+mintermedia2)*c1
       mintermedia4 <- (WM*c2-mintermedia3)^2
       Ow4 <- sum(ConfM*mintermedia4)
       K <- (Ow1-Ow2)/c2
       SdK <- sqrt((Ow4-(Ow1*Ow2-2*Ow2+Ow1)^2)/(SumaMatriz*(c2^4)))
       CV <- K/SdK
     return(list(WeightMatrix=WM, K=K, SdK=SdK, CV=CV))
     },


# Tau ---------------------------------------------------------------------



      #' @description Public method that calculates the Tau index and
      #' its variance. Its value indicates how much the classification has
      #' improved compared to a random classification of the N elements into
      #' M groups. The method also offers the
      #' variance and confidence interval.
      #' The reference \insertCite{ma1995Tau;textual}{ConfMatrix} is followed
      #' for the computations.
      #'
      #' \deqn{
      #' Tau = \dfrac{OverallAcc-PrAgCoef}{1-PrAgCoef}
      #' }
      #' \deqn{
      #' PrAgCoef=\dfrac{1}{M}
      #' }
      #'
      #' \deqn{
      #' \sigma^2_{Tau}=\dfrac{OverallAcc \cdot (1-OverallAcc)}
      #' {N_{Total} \cdot (1-PrAgCoef)^2}
      #' }
      #'
      #' Where:
      #' \enumerate{
      #'   \item \eqn{OverallAcc}: overall accuracy.
      #'   \item \eqn{PrAgCoef}: a priori random agreement coefficient.
      #' }
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #' @return A list of real values containing the Tau index,
      #' its variance and confidence interval.
      #'
      #' @examples
      #' A<-matrix(c(238051,7,132,0,0,24,9,2,189,1,4086,188,0,4,16,45,1,0,939,5082,
      #' 51817,0,34,500,1867,325,17,0,0,5,11148,1618,78,0,0,0,0,48,4,834,2853,340,
      #' 32,0,197,5,151,119,135,726,6774,75,1,553,0,105,601,110,174,155,8257,8,0,
      #' 29,36,280,0,0,6,5,2993,0,115,2,0,4,124,595,0,0,4374),nrow=9,ncol=9)
      #' p<-ConfMatrix$new(A,Source="Muñoz 2016")
      #' p$Tau()
      #'
      #' @aliases NULL

     Tau = function(a=NULL){
        Ca<-1/nrow(self$Values)
        Tau <- ((self$OverallAcc()[[1]]-Ca)/(1-Ca))
        VarTau <- ((self$OverallAcc()[[1]]*(1-self$OverallAcc()[[1]]))/(sum(self$Values)*(1-Ca)^2))
        ConfInt <- private$ConfInt(Tau,VarTau,a)
     return(list(Tau=Tau,VarTau=VarTau,
                 Conf_Int=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup)))
     },



      #' @description Public method that calculates the general Tau
      #' concordance index (weighted) and its standard deviation.
      #'
      #' Be \eqn{p_{ij}=\dfrac{x_{ij}}{\sum^M_{i,j} x_{ij}}} for each element \eqn{i,j} for the matrix
      #' and \eqn{0 \leq w_{ij} \leq 1} for \eqn{i \neq j} and \eqn{w_{ii}=1} for \eqn{i = j}.
      #' If the elements of the weight are greater than 1, their value must be given as a percentage.
      #' The following real values are defined:
      #'
      #' \deqn{
      #' \theta_1=\sum_{i}^{M} p_{ii}
      #' }
      #'
      #' \deqn{
      #' \theta_2=\sum^M_{i=1} w_{ij}p_{i+}
      #' }
      #'
      #'
      #' \deqn{
      #' \theta_3=\sum^M_{i=1} \left( p_{ii} (w_{ij}+p_{+j}) \right)
      #' }
      #'
      #'
      #' \deqn{
      #' \theta_4=\sum^M_{i,j=1} p_{ij} m_{ij}
      #' }
      #'
      #' where \eqn{m_{ij}} are the elements of a matrix, which are given by \eqn{(w_{ij}+p_{+j})^2}
      #'
      #' Therefore,
      #'
      #' \deqn{
      #' Tau_W=\dfrac{\theta_1-\theta_2}{1-\theta_2}
      #' }
      #'
      #' \deqn{
      #' \sigma^2_{Tau_W}=\dfrac{1}{N_{Total}} \left( \dfrac{\theta_1 (1-\theta_1)}{(1-\theta_2)^2}
      #'  + 2 \dfrac{1-\theta_1}{(1-\theta_2)^3} (2 \theta_1 \theta_2-\theta_3) +
      #'  \dfrac{(1-\theta_1)^2}{(1-\theta_2)^4} (\theta_4 - 4 \theta_2^2) \right)
      #' }
      #'
      #' The statistic is given by
      #'
      #'  \deqn{
      #' Z=\dfrac{Tau_W}{\sqrt{\sigma^2_{Tau_W}}}
      #' }
      #'
      #' Where:
      #' \enumerate{
      #'   \item \eqn{\theta_1, \theta_2, \theta_3, \theta_4}: real values.
      #'   \item \eqn{Z}: the test statistic.
      #' }
      #'
      #' @param WV \verb{
      #' Weights vector (as matrix)
      #' }
      #' @return  A list with the weighted Tau index, the weight matrix,
      #' its standard deviation and its statistics.
      #' @examples
      #' A<-matrix(c(65,6,0,4,4,81,11,7,22,5,85,3,24,8,19,90),nrow=4,ncol=4)
      #' WV <-matrix(c(0.4, 0.1, 0.4, 0.1),ncol=4)
      #' p<-ConfMatrix$new(A,Source="Congalton and Green 2008")
      #' p$DetailWTau(WV)
      #'
      #' @aliases NULL

     DetailWTau = function(WV){
       nc <- nrow(self$Values)
       SumaMatriz <-sum(self$Values)
       ConfM<- self$Values/SumaMatriz
       pcol <- apply(ConfM,2,sum)
       prow<- apply(ConfM,1,sum)
       O1 <- sum(diag(ConfM)  )
       O2 <- sum(WV*pcol)
       O3 <- sum(diag(ConfM)*(WV+pcol))
       mintermedia1<- matrix(rep(pcol, nc), nrow =nc, ncol=nc, byrow=FALSE)
       mintermedia2<- matrix(rep(WV, nc), nrow =nc, ncol=nc, byrow=TRUE)
       mintermedia3 <-(mintermedia1+mintermedia2)^2
       O4 <- sum(ConfM*mintermedia3)
       t1<- (1-O1) 
       t2<- (1-O2) 
       t3<- O1*t1/(t2^2)
       t4<- 2*t1*(2*O1*O2-O3)/(t2^3)
       t5<- (t1^2)*(O4-4*O2^2)/(t2^4)
       Tau <- (O1-O2)/t2
       SdT <- sqrt((t3+t4+t5)/SumaMatriz)
       CV<- Tau/SdT
     return(list(Tau=Tau, WeightsVector=WV, SdT=SdT, CV=CV))
     },


# Entropy -----------------------------------------------------------------



      #' @description Public method for calculating product entropy,which
      #' refers to the lack of orden and predictability that the product
      #' presents. The method also offers the variance and confidence
      #' interval. The reference \insertCite{finn1993;textual}{ConfMatrix} is
      #' followed for the calculations.
      #'  \deqn{
      #' Ent=\sum^M_{i,j=1} \left(\dfrac{x_{ij}}{\sum^M_{i,j=1} x_{ij}}
      #'  \cdot \log \left(\dfrac{x_{ij}}{\dfrac{ x_{i+}
      #'  \cdot  x_{+j}}{\sum^M_{i,j=1} x_{ij}}} \right) \right)
      #' }
      #' \deqn{
      #' \sigma^2_{Ent}=\dfrac{Ent \cdot (1-Ent)}{N_{Total}}
      #' }
      #'
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #'
      #' @param v \verb{
      #' Base of the logarithm, where} \eqn{v \in \mathbb{R}^{+}-\{1\}}. \verb{By default v=10(units Hartleys),
      #' v=2(units bits), v=e(units nats).}
      #' @return A list of real values containing the entropy, its variance
      #' and confidence interval.
      #' @examples
      #' A<-matrix(c(35,4,12,2,14,11,9,5,11,3,38,12,1,0,4,2),nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A,Source="Finn 1993")
      #' p$Ent(v=2)
      #'
      #' @aliases NULL

     Ent = function(a=NULL,v=NULL){
      if(!is.null(v)){
        v<-v
      }else{v<-10}

       pbi<-private$sumfil(self$Values)/sum(self$Values)
       paj<-private$sumcol(self$Values)/sum(self$Values)
       pbiaj<-self$Values/sum(self$Values)
       pbiaj_<-matrix(rep(0),nrow=nrow(self$Values),ncol=ncol(self$Values))
       for (i in 1:nrow(self$Values)) {
         pbiaj_[,i]<-pbiaj[,i]/paj[i]
       }
       pbiaj_pbi<-matrix(rep(0),nrow=nrow(self$Values),ncol=ncol(self$Values))
       for (j in 1:nrow(self$Values)) {
         pbiaj_pbi[j,]<-pbiaj_[j,]/pbi[j]
       }

       l<-log(pbiaj_pbi,base=v)
       s<-pbiaj*l
       s<-s[!is.na(s)]
       Ent<-sum(s)

       VarEnt <- abs((Ent*(1-Ent))/sum(self$Values))
       ConfInt <- private$ConfInt(Ent,VarEnt,a)
     return(list(Entropy=Ent,VarEnt=VarEnt,
                 Conf_Int=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup)))
     },




      #' @description Public method that calculates normalized entropy using
      #' the arithmetic mean of the entropies on the product and the
      #' reference. The method also offers the variance and confidence interval. The reference
      #' \insertCite{strehl2002;textual}{ConfMatrix} is followed for the calculations.
      #' \deqn{
      #' AvNormEnt=\dfrac{2Ent}{Ent_i(A)+Ent_i(B)}
      #' }
      #'  \deqn{
      #' Ent_i(A)=-\sum^M_{j=1} \left( \left(\dfrac{ x_{+j}}
      #' {\sum^M_{i,j=1} x_{ij} }\right) \cdot \log \left(\dfrac{ x_{+j}}
      #' {\sum^M_{i,j=1} x_{ij} }\right) \right)
      #' }
      #' \deqn{
      #' Ent_i(B)=-\sum^M_{i=1}\left( \left(\dfrac{ x_{i+}}
      #' {\sum^M_{i,j=1} x_{ij} }\right) \cdot \log \left(\dfrac{x_{i+}}
      #' {\sum^M_{i,j=1} x_{ij} }\right) \right)
      #' }
      #' \deqn{
      #' \sigma^2_{AvNormEnt}=\dfrac{AvNormEnt \cdot (1-AvNormEnt)}{N_{Total}}
      #' }
      #'
      #' where:
      #'
      #' \enumerate{
      #'   \item \eqn{Ent}: product entropy.
      #'   \item \eqn{Ent_i(A)}: entropy with respect to the classes \emph{i}
      #'   of the product. A is a matrix.
      #'   \item \eqn{Ent_i(B)}: entropy with respect to the class \emph{i} on the reference. B is a matrix.
      #' }
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #' @param v \verb{
      #' Base of the logarithm, where} \eqn{v \in \mathbb{R}^{+}-\{1\}}. \verb{By default v=10(units Hartleys),
      #' v=2(units bits), v=e(units nats).}
      #' @return A list of real values containing the normalized
      #' entropy (arithmetic mean of the entropies on the product
      #' and reference), its variance and confidence interval.
      #' @examples
      #' A<-matrix(c(0,12,0,0,12,0,0,0,0,0,0,12,0,0,12,0),nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A,Source="Liu et al. 2007")
      #' p$AvNormEnt(v=2)
      #'
      #' @aliases NULL

     AvNormEnt = function(a=NULL,v=NULL){
      if(!is.null(v)){
        v<-v
      }else{v<-10}

      Ent_iB <- - sum ((private$sumfil(self$Values)/sum(self$Values)) *
                     (log(private$sumfil(self$Values)/sum(self$Values),base=v)),na.rm=TRUE)
      Ent_iA <- - sum ((private$sumcol(self$Values)/sum(self$Values)) *
                     (log(private$sumcol(self$Values)/sum(self$Values),base=v)),na.rm=TRUE)
        if (Ent_iA + Ent_iB == 0) {
          stop ("/ by 0")
        }else{
          AvNormEnt <- 2 * self$Ent(v=v)[[1]] / (Ent_iA + Ent_iB)
          VarAvNormEnt <- abs((AvNormEnt*(1-AvNormEnt))/sum(self$Values))
          ConfInt <- private$ConfInt(AvNormEnt,VarAvNormEnt,a)
        }

     return(list(AvNormEntrop=AvNormEnt,VarAvNormEntrop=VarAvNormEnt,
                 Conf_Int=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup)))
     },



      #' @description Public method that calculates the normalized entropy
      #' using the geometric mean of the product and reference entropies.
      #' The method also offers the variance and confidence interval.
      #' The reference \insertCite{ghosh2002;textual}{ConfMatrix} is followed
      #' for the calculations.
      #' @description
      #'
      #' \deqn{
      #' GeomAvNormEnt=\dfrac{Ent}{\sqrt{Ent_i(A) \cdot Ent_i(B)}}
      #' }
      #'  \deqn{
      #' Ent_i(A)=-\sum^M_{j=1} \left( \left(\dfrac{ x_{+j}}
      #' {\sum^M_{i,j=1} x_{ij} }\right) \cdot \log \left(\dfrac{ x_{+j}}
      #' {\sum^M_{i,j=1} x_{ij} }\right) \right)
      #' }
      #' \deqn{
      #' Ent_i(B)=-\sum^M_{i=1}\left( \left(\dfrac{ x_{i+}}
      #' {\sum^M_{i,j=1} x_{ij} }\right) \cdot \log \left(\dfrac{ x_{i+}}
      #' {\sum^M_{i,j=1} x_{ij} }\right) \right)
      #' }
      #'\deqn{
      #' \sigma^2_{GeomAvNormEnt}=
      #' \dfrac{GeomAvNormEnt \cdot (1-GeomAvNormEnt)}{N_{Total}}
      #' }
      #' where:
      #'
      #' \enumerate{
      #'   \item \eqn{Ent}: product entropy.
      #'   \item \eqn{Ent_i(A)}: entropy with respect to the classes \emph{i}
      #'   of the product. A is a matrix.
      #'   \item \eqn{Ent_i(B)}: entropy with respect to the class \emph{i} of the reference. B is a matrix.
      #' }
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #' @param v \verb{
      #' Base of the logarithm, where} \eqn{v \in \mathbb{R}^{+}-\{1\}}. \verb{By default v=10(units Hartleys),
      #' v=2(units bits), v=e(units nats).}
      #' @return A list of real values containing the normalized
      #' entropy (geometric mean of the entropies on the product
      #' and reference), its variance and confidence interval.
      #'
      #' @examples
      #' A<-matrix(c(0,12,0,0,12,0,0,0,0,0,0,12,0,0,12,0),nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A,Source="Liu et al. 2007")
      #' p$GeomAvNormEnt(v=2)
      #'
      #' @aliases NULL

     GeomAvNormEnt = function(a=NULL,v=NULL){
      if(!is.null(v)){
        v<-v
      }else{v<-10}
      Ent_iB <- - sum ((private$sumfil(self$Values)/sum(self$Values)) * (log(private$sumfil(self$Values)/sum(self$Values),base=v)),na.rm=TRUE)
      Ent_iA <- - sum ((private$sumcol(self$Values)/sum(self$Values)) * (log(private$sumcol(self$Values)/sum(self$Values),base=v)),na.rm=TRUE)
       if (Ent_iA * Ent_iB == 0) {
        stop ("/ by 0")
       }else{
         GeomAvNormEnt <- self$Ent(v=v)[[1]] / sqrt(Ent_iA * Ent_iB)
         VarGeomAvNormEnt <- abs((GeomAvNormEnt*(1-GeomAvNormEnt))/sum(self$Values))
         ConfInt<-private$ConfInt(GeomAvNormEnt,VarGeomAvNormEnt,a)
         }
     return(list(GeomAvNormEntrop=GeomAvNormEnt,
                 VarGeoAvNormEntrop=VarGeomAvNormEnt,
                 Conf_Int=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup)))
     },



      #' @description Public method that provides normalized entropy using
      #' the arithmetic mean of the maximum entropies of the product and
      #' reference. The method also offers the variance and confidence interval.
      #' The reference \insertCite{strehl2002relationship;textual}{ConfMatrix} is
      #' followed for the calculations.
      #' @description
      #' \deqn{
      #' AvMaxNormEnt=\dfrac{2 Ent}{max(Ent_i(A))+max(Ent_i(B))}=
      #' \dfrac{Ent}{\log M}
      #' }
      #'  \deqn{
      #' Ent_i(A)=-\sum^M_{j=1} \left( \left(\dfrac{ x_{+j}}
      #' {\sum^M_{i,j=1} x_{ij} }\right) \cdot \log \left(\dfrac{ x_{+j}}
      #' {\sum^M_{i,j=1} x_{ij} }\right) \right)
      #' }
      #' \deqn{
      #' Ent_i(B)=-\sum^M_{i=1}\left( \left(\dfrac{ x_{i+}}
      #' {\sum^M_{i,j=1} x_{ij} }\right) \cdot \log \left(\dfrac{ x_{i+}}
      #' {\sum^M_{i,j=1} x_{ij} }\right) \right)
      #' }
      #' \deqn{
      #' \sigma^2_{AvMaxNormEnt}=
      #' \dfrac{AvMaxNormEnt \cdot (1-AvMaxNormEnt)}{N_{Total}}
      #' }
      #' where:
      #'
      #' \enumerate{
      #'   \item \eqn{Ent}: product entropy.
      #'   \item \eqn{Ent_i(A)}: entropy with respect to the classes \emph{i}
      #'   of the product. A is a matrix.
      #'   \item \eqn{Ent_i(B)}: entropy with respect to the class \emph{i} on the reference. B is a matrix.
      #' }
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #' @param v \verb{
      #' Base of the logarithm, where} \eqn{v \in \mathbb{R}^{+}-\{1\}}. \verb{By default v=10(units Hartleys),
      #' v=2(units bits), v=e(units nats).}
      #' @return A list of real values containing the normalized entropy
      #' (arithmetic mean of the maximum entropies of the product and of
      #' reference), its variance, and its confidence interval.
      #'
      #' @examples
      #' A<-matrix(c(8,0,0,0,0,16,0,0,0,0,8,0,0,0,0,16),nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A,Source="Liu et al. 2007")
      #' p$AvMaxNormEnt(v=2)
      #'
      #' @aliases NULL

     AvMaxNormEnt = function(a=NULL,v=NULL){
       if(!is.null(v)){
         v<-v
       }else{v<-10}

        AvMaxNormEnt <- self$Ent(v=v)[[1]] / log(sqrt(length(self$Values)),base=v)
        VarAvMaxNormEnt <- abs((AvMaxNormEnt*(1-AvMaxNormEnt))/sum(self$Values))
        ConfInt <- private$ConfInt(AvMaxNormEnt,VarAvMaxNormEnt,a)

     return (list(AvMaxNormEnt=AvMaxNormEnt,Var=VarAvMaxNormEnt,
                  Conf_Int=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup)))
     },



     #' @description Public method that calculates relative change of
     #' entropy for a given class \eqn{i} of the product. The method also
     #' offers the variance and confidence interval.
     #' The reference \insertCite{finn1993;textual}{ConfMatrix} is followed for
     #' the calculations.
     #' @description
     #' \deqn{
     #' EntUser_i= \dfrac{Ent_i(A)-Ent_i(A|b_i)}{Ent_i(A)}
     #' }
     #'  \deqn{
     #' Ent_i(A)=-\sum^M_{j=1} \left( \left(\dfrac{ x_{+j}}
     #' {\sum^M_{i,j=1} x_{ij} }\right) \cdot \log \left(\dfrac{ x_{+j}}
     #' {\sum^M_{i,j=1} x_{ij} }\right) \right)
     #' }
     #' \deqn{
     #' Ent_i(A|b_i)=-\sum^M_{j=1} \left( \left(\dfrac{ x_{ij}}
     #' { x_{i+} }\right) \cdot \log \left(\dfrac{x_{ij}}
     #' { x_{i+}}\right) \right)
     #' }
     #' \deqn{
     #' \sigma^2_{EntUser_i}= \dfrac{EntUser_i \cdot (1-EntUser_i)}{N_{Total}}
     #' }
     #'
     #' where:
     #'
     #' \enumerate{
     #'   \item \eqn{Ent_i(A)}: entropy with respect to the classes \emph{i}
     #'   of the product. A is a matrix.
     #'   \item \eqn{Ent_i(A|b_i)}: Producer entropy knowing that the
     #'   location corresponding to reference B is in class \eqn{b_i}.
     #'   B is a matrix.
     #' }
     #' @param i \verb{
     #' Class to evaluate, where} \eqn{i \in \mathbb{Z}-\{0\}}.
     #'
     #' @param a \verb{
     #' Significance level. By default 0.05.
     #' }
     #' @param v \verb{
     #' Base of the logarithm, where} \eqn{v \in \mathbb{R}^{+}-\{1\}}. \verb{By default v=10(units Hartleys),
     #' v=2(units bits), v=e(units nats).}
     #'
     #' @return  A list of real values containing the relative change of entropy
     #' for given class i, its variance, its confidence interval, producer's
     #' entropy, and producer's entropy knowing that the
     #' location corresponding to reference B is in class \eqn{b_i}.
     #'
     #' @examples
     #' A<-matrix(c(35,4,12,2,14,11,9,5,11,3,38,12,1,0,4,2),nrow=4,ncol=4)
     #' p<-ConfMatrix$new(A,Source="Finn 1993")
     #' p$EntUser_i(1,v=2)
     #'
     #' @aliases NULL

    EntUser_i = function(i,a=NULL,v=NULL){
      if(!is.null(v)){
       v<-v
      }else{v<-10}

      Ent_iA <- - sum ((private$sumcol(self$Values)/sum(self$Values)) *
                  (log(private$sumcol(self$Values)/sum(self$Values),base=v)),na.rm=TRUE)
      Ent_iAbi <- - sum ((self$Values[i,] / private$sumfil(self$Values)[i]) *
                    log(self$Values[i,] / private$sumfil(self$Values)[i],base=v),na.rm=TRUE)

      if (Ent_iA == 0){
        stop("/by 0")
      }else {
        EntUser_i <- (Ent_iA - Ent_iAbi) / Ent_iA
        VarEntUser_i <- abs((EntUser_i*(1-EntUser_i))/sum(self$Values))
        ConfInt <- private$ConfInt(EntUser_i,VarEntUser_i,a)
      }
      return(list(EntUser_i=EntUser_i,VarEntUser_i=VarEntUser_i,
                  Conf_Int=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup),
                  Entrop_iA=Ent_iA,Entrop_iAbi=Ent_iAbi))
    },



      #' @description Public method that calculates normalized entropy
      #' of the product. The method also offers the variance and
      #' confidence interval. The reference
      #' \insertCite{finn1993;textual}{ConfMatrix} is followed for the calculations.
      #' @description
      #'
      #' \deqn{
      #' NormEntUser=\dfrac{Ent}{Ent_i(B)}
      #' }
      #' \deqn{
      #' Ent_i(B)=-\sum^M_{i=1} \left( \left(  \dfrac {x_{i+}}
      #' {\sum^M_{i,j=1} x_{ij} }\right) \cdot \log \left( \dfrac{ x_{i+}}
      #' {\sum^M_{i,j=1} x_{ij} }\right) \right)
      #' }
      #'\deqn{
      #' \sigma^2_{NormEntUser}=\dfrac{NormEntUser \cdot (1-NormEntUser)}{N_{Total}}
      #' }
      #' where:
      #'
      #' \enumerate{
      #'   \item \eqn{Ent}: product entropy.
      #'   \item \eqn{Ent_i(B)}: entropy with respect to the class \emph{i} on the reference. B is a matrix.
      #' }
      #'
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #'
      #' @param v \verb{
      #' Base of the logarithm, where} \eqn{v \in \mathbb{R}^{+}-\{1\}}. \verb{By default v=10(units Hartleys),
      #' v=2(units bits), v=e(units nats).}
      #' @return A list of real values containing with normalized entropy
      #' of the product class i, conditioned to reference data, its variance
      #' and confidence interval.
      #'
      #' @examples
      #' A<-matrix(c(35,4,12,2,14,11,9,5,11,3,38,12,1,0,4,2),nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A,Source="Finn 1993")
      #' p$NormEntUser(v=2)
      #'
      #' @aliases NULL

     NormEntUser = function(a=NULL,v=NULL){
      if(!is.null(v)){
        v<-v
      }else{v<-10}

       Ent_iB <- - sum ((private$sumfil(self$Values)/sum(self$Values)) *
                      (log(private$sumfil(self$Values)/sum(self$Values),base=v)),na.rm=TRUE)

       if(Ent_iB == 0){
        stop("/ by 0")
       }
       NormEntUser <- self$Ent(v=v)[[1]]/Ent_iB
       VarNormEntUser <- abs((NormEntUser*(1-NormEntUser))/sum(self$Values))
       ConfInt <- private$ConfInt(NormEntUser,VarNormEntUser,a)
     return(list(NormEntropUser=NormEntUser,
                 VarNormEntropUser=VarNormEntUser,
                 Conf_Int=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup)))
     },




     #' @description Public method that calculates relative change of
     #' entropy for a given a class \eqn{i} of the reference from the producer's
     #' perspective. The method also offers the variance and confidence interval.
     #' The reference \insertCite{stehman1997;textual}{ConfMatrix} is followed for
     #' the calculations.
     #'
     #' \deqn{
     #' EntProd_i= \dfrac{Ent_i(B)-Ent_i(B|a_j)}{Ent_i(B)}
     #' }
     #'
     #'  \deqn{
     #' Ent_i(B)=-\sum^M_{i=1} \left( \left(\dfrac{ x_{i+}}
     #' {\sum^M_{i,j=1} x_{ij} }\right) \cdot \log \left(\dfrac{ x_{i+}}
     #' {\sum^M_{i,j=1} x_{ij} }\right) \right)
     #' }
     #' \deqn{
     #' Ent_i(B|a_j)=-\sum^M_{i=1}\left( \left(\dfrac{ x_{ij}}
     #' { x_{+j} }\right) \cdot \log \left(\dfrac{x_{ij}}
     #' { x_{+j}}\right) \right)
     #' }
     #'\deqn{
     #' \sigma^2_{EntProd_i}= \dfrac{EntProd_i \cdot (1-EntProd_i)}{N_{Total}}
     #' }
     #' where:
     #'
     #' \enumerate{
     #'   \item \eqn{Ent_i(B)}: entropy with respect to the class \emph{i} on the reference. B is a matrix.
     #'   \item \eqn{Ent_i(B|a_j)}: Entropy of reference B knowing that the
     #'   location of product A is in the class \eqn{a_j}.
     #' }
     #'
     #' @param i \verb{
     #' Class to evaluate, where} \eqn{i \in \mathbb{Z}-\{0\}}.
     #'
     #' @param a \verb{
     #' Significance level. By default 0.05.
     #' }
     #' @param v \verb{
     #' Base of the logarithm, where} \eqn{v \in \mathbb{R}^{+}-\{1\}}. \verb{ By default v=10(units Hartleys),
     #' v=2(units bits), v=e(units nats).
     #' }
     #' @return
     #' A list of real values containing the relative change of entropy
     #' for given class i, its variance, its confidence interval, entropy with
     #' respect to reference classes, and entropy with respect to reference
     #' classes knowing that the location corresponding to A is in class \eqn{a_j}.
     #' @examples
     #' A<-matrix(c(35,4,12,2,14,11,9,5,11,3,38,12,1,0,4,2),nrow=4,ncol=4)
     #' p<-ConfMatrix$new(A,Source="Finn 1993")
     #' p$EntProd_i(3,v=2)
     #'
     #' @aliases NULL

    EntProd_i = function(i,a=NULL,v=NULL){
     if(!is.null(v)){
      v<-v
     }else{v<-10}
    Ent_iB <- - sum ((private$sumfil(self$Values)/sum(self$Values)) *
                  (log(private$sumfil(self$Values)/sum(self$Values),base = v)),na.rm=TRUE)
    Ent_iBaj <- - sum ((self$Values[,i] / private$sumcol(self$Values)[i]) *
                    log(self$Values[,i] / private$sumcol(self$Values)[i],base=v),na.rm=TRUE)

      if (Ent_iB == 0){
        stop("/by 0")
      }else {
        EntProd_i <- (Ent_iB - Ent_iBaj) / Ent_iB
        VarEntProd_i <- abs((EntProd_i*(1-EntProd_i))/sum(self$Values))
        ConfInt <- private$ConfInt(EntProd_i,VarEntProd_i,a)
        }
    return(list(EntProd_i=EntProd_i,VarEntProd_i=VarEntProd_i,
                Conf_Int=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup),
                Entrop_iB=Ent_iB,Entrop_iBaj=Ent_iBaj))
    },



      #' @description Public method that calculates normalized entropy of
      #' the reference from the producer's perspective. The method also offers the variance and confidence
      #' interval. The reference \insertCite{finn1993;textual}{ConfMatrix} is
      #' followed for the calculations.
      #' @description
      #' \deqn{
      #' NormEntProd=\dfrac{Ent}{Ent_i(A)}
      #' }
      #'  \deqn{
      #' Ent_i(A)=-\sum^M_{j=1}\left( \left(\dfrac{ x_{+j}}
      #' {\sum^M_{i,j=1} x_{ij} }\right) \cdot \log \left( \dfrac{x_{+j}}
      #' {\sum^M_{i,j=1} x_{ij} }\right) \right)
      #' }
      #'\deqn{
      #' \sigma^2_{NormEntProd}=\dfrac{NormEntProd \cdot (1-NormEntProd)}{N_{Total}}
      #' }
      #' where:
      #'
      #' \enumerate{
      #'   \item \eqn{Ent}: product entropy.
      #'   \item \eqn{Ent_i(A)}: entropy with respect to the classes \emph{i}
      #'   of the product. A is a matrix.
      #' }
      #'
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #'
      #' @param v \verb{
      #' Base of the logarithm, where} \eqn{v \in \mathbb{R}^{+}-\{1\}}.\verb{ By default v=10(units Hartleys),
      #' v=2(units bits), v=e(units nats).}
      #' @return A list of real values containing the normalized entropy
      #' of the reference class \eqn{i} from the producer's perspective, its variance
      #' and confidence interval.
      #'
      #' @examples
      #' A<-matrix(c(35,4,12,2,14,11,9,5,11,3,38,12,1,0,4,2),nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A,Source="Finn 1993")
      #' p$NormEntProd(v=2)
      #'
      #' @aliases NULL

     NormEntProd = function(a=NULL,v=NULL){
      if(!is.null(v)){
        v<-v
      }else{v<-10}

      Ent_iA <- - sum ((private$sumcol(self$Values)/sum(self$Values)) *
                     (log(private$sumcol(self$Values)/sum(self$Values),base=v)),na.rm=TRUE)
       if (Ent_iA == 0){
        stop ("/ by 0")
       }else{
         NormEntProd <- self$Ent(v=v)[[1]]/Ent_iA
         VarNormEntProd <- abs((NormEntProd*(1-NormEntProd))/sum(self$Values))
         ConfInt <- private$ConfInt(NormEntProd,VarNormEntProd,a)
         }
     return(list(NormEntropProd=NormEntProd,
            VarNormEntropProd=VarNormEntProd,
            Conf_Int=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup)))
     },


# Others index -------------------------------------------------------------------------


      #' @description Public method that provides the Classification
      #' Success Index (CSI) which considers all classes and gives an
      #' overall estimation of classification effectiveness.
      #' The method also offers the variance and confidence interval.
      #' The references \insertCite{koukoulas2001;textual}{ConfMatrix} and
      #' \insertCite{turk2002;textual}{ConfMatrix} are followed for the calculations.
      #'  \deqn{
      #' Sucess=1-(1-AvUserAcc+1-AvProdAcc)=AvUserAcc+AvProdAcc-1
      #' }
      #'  \deqn{
      #' \sigma^2_{Sucess}=\dfrac{Sucess \cdot (1-Sucess)}{N_{Total}}
      #' }
      #'
      #' where:
      #'
      #' \enumerate{
      #'   \item \eqn{AvUserAcc}: average accuracy from user's perspective.
      #'   \item \eqn{AvProdAcc}: average accuracy from producer's perspective.
      #' }
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #' @return A list of real values containing the ICSI, its variance
      #' and its confidence interval.
      #' @examples
      #' A<-matrix(c(0.3,0.02,0.01,0.12,0.19,0.03,0.02,0.01,0.3),nrow=3,ncol=3)
      #' p<-ConfMatrix$new(A,Source="Labatut and Cherifi 2011")
      #' p$Sucess()
      #'
      #' @aliases NULL

     Sucess = function(a=NULL){
      Sucess <- self$AvUserAcc()[[1]] + self$AvProdAcc()[[1]] - 1
       VarSucess <- abs((Sucess*(1-Sucess))/sum(self$Values))
       ConfInt <- private$ConfInt(Sucess,VarSucess,a)
     return(list(Sucess=Sucess,VarSucess=VarSucess,
                 Conf_Int=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup)))
     },



      #' @description Public method that provides the Individual
      #' Classification Success Index (ICSI) which considers  the classification
      #' effectiveness for one particular class of interest.
      #' The method also offers the variance and confidence interval.
      #' The references \insertCite{koukoulas2001;textual}{ConfMatrix} and \insertCite{turk2002;textual}{ConfMatrix}
      #' are followed for the calculations.
      #'  \deqn{
      #' Sucess_i=1-(1-UserAcc_i+1-ProdAcc_i)=UserAcc_i+ProdAcc_i-1
      #' }
      #'
      #' \deqn{
      #' \sigma^2_{Sucess_i}=\dfrac{Sucess_i \cdot (1-Sucess_i)}{N_{ij}}
      #' }
      #'
      #' where:
      #'
      #' \enumerate{
      #'   \item \eqn{UserAcc_i}: user accuracy index for class i.
      #'   \item \eqn{ProdAcc_i}: producer accuracy index for class i.
      #' }
      #'
      #' @param i \verb{
      #' Class to evaluate, where} \eqn{i \in \mathbb{Z}-\{0\}}.
      #'
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #' @return A list of real values containing the ICSI,
      #' its variance and its confidence interval.
      #'
      #' @examples
      #' A<-matrix(c(0.3,0.02,0.01,0.12,0.19,0.03,0.02,0.01,0.3),nrow=3,ncol=3)
      #' p<-ConfMatrix$new(A,Source="Labatut and Cherifi 2011")
      #' p$Sucess_i(2)
      #'
      #' @aliases NULL

     Sucess_i = function(i,a=NULL){
      Sucess_i <- self$UserAcc_i(i)[[1]] + self$ProdAcc_i(i)[[1]] - 1
      VarSucess_i <- abs((Sucess_i*(1-Sucess_i))/(private$sumcol(self$Values)[i]+private$sumfil(self$Values)[i]-self$Values[i,i]))
      ConfInt <- private$ConfInt(Sucess_i,VarSucess_i,a)
     return (list(Sucess_i=Sucess_i,VarSucess_i=VarSucess_i,
                  Conf_Int=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup)))
     },



      #' @description Public method that provides the average value of the
      #' Hellden mean precision index. Denoted by the probability that a
      #' randomly chosen position or element assigned to a specific class on
      #' the product has a correspondence of the same class in the homologous
      #' position or element in the reference, and that a randomly chosen point
      #' or element assigned to a specific class on the reference has a
      #' correspondence of the same class in the homologous position or
      #' element in the product. The method also offers the
      #' variance and confidence interval.
      #' The reference \insertCite{liu2007;textual}{ConfMatrix} is followed for
      #' the calculations.
      #' @description
      #'  \deqn{
      #' AvHellAcc=\dfrac{1}{M} 2 \sum^M_{i=1} \dfrac{ x_{ii}}
      #' { x_{+i} + x_{i+}}
      #' }
      #'  \deqn{
      #' \sigma^2_{AvHellAcc}=\dfrac{AvHellAcc \cdot (1-AvHellAcc)}{N_{Total}}
      #' }
      #'
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #' @return A list of real values containing the average of
      #' Hellden's mean accuracy index, its variance and
      #' confidence interval.
      #'
      #' @examples
      #' A<-matrix(c(65,6,0,4,4,81,11,7,22,5,85,3,24,8,19,90),nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A,Source="Congalton and Green 2008")
      #' p$AvHellAcc()
      #'
      #' @aliases NULL

     AvHellAcc = function(a=NULL){
      AvHellAcc <- 1/sqrt(length(self$Values)) *
        sum ((2*diag(self$Values)) / (private$sumfil(self$Values) + private$sumcol(self$Values)))
      VarAvHellAcc <- abs((AvHellAcc*(1-AvHellAcc))/sum(self$Values))
      ConfInt <- private$ConfInt(AvHellAcc,VarAvHellAcc,a)
     return(list(AvHellAcc=AvHellAcc,VarAvHellAcc=VarAvHellAcc,
                 Conf_Int=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup)))
     },



      #' @description Public method that provides the Hellden’ average
      #' accuracy for the specified class. The method also offers the variance and
      #' confidence interval. The references
      #' \insertCite{hellden1980;textual}{ConfMatrix} and \insertCite{rosenfield1986;textual}{ConfMatrix} are
      #' followed for the calculations.
      #' @description
      #'  \deqn{
      #' AvHellAcc_i=\dfrac{2}{\dfrac{1}{UserAcc_i}+\dfrac{1}{ProdAcc_i}}=
      #' \dfrac{2 UserAcc_i \cdot ProdAcc_i}{UserAcc_i + ProdAcc_i}
      #' }
      #' \deqn{
      #' \sigma^2_{AvHellAcc_i}=\dfrac{AvHellAcc_i \cdot (1-AvHellAcc_i)}{N_{ij}}
      #' }
      #' where:
      #'
      #' \enumerate{
      #'   \item \eqn{UserAcc_i}: user accuracy index for class i.
      #'   \item \eqn{ProdAcc_i}: producer accuracy index for class i.
      #' }
      #' @param i \verb{
      #' Class to evaluate, where} \eqn{i \in \mathbb{Z}-\{0\}}.
      #'
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #' @return A list of real values containing the Hellden’s mean
      #' accuracy, its variance and its confidence interval.
      #'
      #' @examples
      #' A <- matrix(c(148,1,8,2,0,0,50,15,3,0,1,6,39,7,1,1,0,6,25,1,1,0,0,1,6),nrow=5,
      #' ncol=5)
      #' p<-ConfMatrix$new(A,Source="Rosenfield and Fitzpatrick 1986")
      #' p$AvHellAcc_i(2)
      #'
      #' @aliases NULL

     AvHellAcc_i = function(i,a=NULL){
       if (self$UserAcc_i(i)[[1]] == 0 || self$ProdAcc_i(i)[[1]] == 0) {
        stop ("/ by 0")
       }else{
          AvHellAcc_i <- 2 / (1/self$UserAcc_i(i)[[1]] +
                                1/self$ProdAcc_i(i)[[1]])
         VarAvHellAcc_i <- abs((AvHellAcc_i*(1-AvHellAcc_i))/
                                    (private$sumcol(self$Values)[i]+private$sumfil(self$Values)[i]-self$Values[i,i]))
         ConfInt <- private$ConfInt(AvHellAcc_i,VarAvHellAcc_i,a)
         }

     return(list(AvHellAcc_i=AvHellAcc_i,
                 VarAvHellAcc_i=VarAvHellAcc_i,
                 Conf_Int=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup)))
     },



      #' @description Public method that provides the average of the Short's
      #' mapping accuracy index. It is determined as the quotient between the
      #' well-classified elements (value on the diagonal) and the subtraction
      #' of that same value on the errors of omission and commission (rest of
      #' values in the column and row) corresponding to each class. The method
      #' also offers the variance and confidence interval. The
      #' reference \insertCite{liu2007;textual}{ConfMatrix} is followed for
      #' the calculations.
      #' @description
      #'  \deqn{
      #' AvShortAcc=\dfrac{1}{M} \sum^M_{i=1} \dfrac{x_{ii}}
      #' { \overline{x}_{+ i}+ \overline{x}_{i +}-x_{ii}}
      #' }
      #'\deqn{
      #' \sigma^2_{AvShortAcc}=\dfrac{AvShortAcc \cdot (1-AvShortAcc)}{N_{Total}}
      #' }
      #'
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #' @return A list of real values containing the average of
      #' Short's mapping accuracy index, its variance and
      #' confidence interval.
      #'
      #' @examples
      #' A<-matrix(c(65,6,0,4,4,81,11,7,22,5,85,3,24,8,19,90),nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A,Source="Congalton and Green 2008")
      #' p$AvShortAcc()
      #'
      #' @aliases NULL

     AvShortAcc = function(a=NULL){
      sum1 <- private$sumfil(self$Values)+private$sumcol(self$Values)- diag(self$Values)
       for (i in 1:length(sum1)) {
          if (sum1[i] == 0) {
           stop ("/ by 0")
          }
       }
      AvShortAcc = 1/sqrt(length(self$Values)) *
        sum (diag(self$Values) / (private$sumfil(self$Values) + private$sumcol(self$Values) - diag(self$Values)))
      VarAvShortAcc=abs((AvShortAcc*(1-AvShortAcc))/sum(self$Values))
      ConfInt <- private$ConfInt(AvShortAcc,VarAvShortAcc,a)

     return(list(AvShortAcc=AvShortAcc,VarAvShortAcc=VarAvShortAcc,
            Conf_Int=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup)))
     },



      #' @description Public method that provides the Short's mapping
      #' accuracy for each class. The method also offers the
      #' variance and confidence interval. The references
      #' \insertCite{rosenfield1986;textual}{ConfMatrix} and \insertCite{short1982;textual}{ConfMatrix}
      #' are followed for the calculations.
      #'  \deqn{
      #' ShortAcc_i=\dfrac{x_{ii}}{ \overline{x}_{+ i}+
      #'  \overline{x}_{i +}-x_{ii}}
      #' }
      #'\deqn{
      #' \sigma^2_{ShortAcc_i}=\dfrac{ShortAcc_i \cdot (1-ShortAcc_i)}{N_{ij}}
      #' }
      #'
      #' @param i \verb{
      #' Class to evaluate, where} \eqn{i \in \mathbb{Z}-\{0\}}.
      #'
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #' @return A list of real values containing the Short's
      #' mapping accuracy, its variance and its confidence interval.
      #'
      #' @examples
      #' A <- matrix(c(148,1,8,2,0,0,50,15,3,0,1,6,39,7,1,1,0,6,25,1,1,0,0,1,6),nrow=5,
      #' ncol=5)
      #' p<-ConfMatrix$new(A,Source="Rosenfield and Fitzpatrick-Lins 1986")
      #' p$ShortAcc_i(2)
      #'
      #' @aliases NULL

     ShortAcc_i = function(i,a=NULL){
      if (private$sumfil(self$Values)[i] + private$sumcol(self$Values)[i] - self$Values[i,i] == 0) {
      stop ("/ by 0")
      }else{
        ShortAcc_i = self$Values[i,i] / (private$sumfil(self$Values)[i] + private$sumcol(self$Values)[i]
                     - self$Values[i,i])
        VarShortAcc_i=abs((ShortAcc_i*(1-ShortAcc_i))/
                      (private$sumcol(self$Values)[i]+private$sumfil(self$Values)[i]-self$Values[i,i]))
        ConfInt <- private$ConfInt(ShortAcc_i,VarShortAcc_i,a)
        }

     return(list(ShortAcc_i=ShortAcc_i,VarShortAcc_i=VarShortAcc_i,
                 Conf_Int=c(ConfInt$ConfInt_inf,ConfInt$ConfInt_sup)))
     },



      #' @description Public method that calculates the Ground Truth index,
      #' its variance and confidence interval.The reference \insertCite{turk1979gt;textual}{ConfMatrix}
      #' is followed for the computations.
      #'
      #' To calculate \eqn{R} we begin the following iterative process:
      #'
      #' Be \eqn{U_j^{(0)}=f_j^0} with \eqn{f_j^0=\dfrac{\overline{x}_{i+}}{\sum_{i=1}^M \overline{x}_{i+}}}
      #' and \eqn{f_i^0=\dfrac{\overline{x}_{+i}}{\sum_{i=1}^M \overline{x}_{+i}}}
      #'
      #'
      #'  Where \eqn{2m} with \eqn{m=1,2,\cdots}
      #'  \deqn{
      #'  V_{i,2m-1}=\dfrac{f_i^0}{U_{+,2m-2}-U_{i,2m-2}}
      #'  }
      #'  where \eqn{U_{+,2m}=\sum_{i=1}^M U_{j,2m} }
      #'  and when \eqn{2m+1} with \eqn{m=1,2,\cdots}
      #'  \deqn{
      #'  U_{j,2m}=\dfrac{f_j^0}{V_{+,2m-1}-V_{i,2m-1}}
      #'  }
      #'  where \eqn{V_{+,2m-1}=\sum_{i=1}^M V_{i,2m-1} }
      #'
      #'  The iterative steps continue for \eqn{m=1, 2,\cdots} until
      #'  the accuracy stabilizes thus taking the V term.
      #'  Where
      #'  \deqn{
      #' R=\dfrac{V}{\sum_{i=1}^{M} V_i}
      #' }
      #'
      #' \deqn{
      #' ProdAcc=\dfrac{x_{ii}}{\sum_{j=1}^M x_{+j}}
      #' }
      #' \deqn{
      #' GroundTruth = \dfrac{ProdAcc-R}{1-R}
      #' }
      #'
      #' \deqn{
      #' \sigma^2_{GroundTruth}=\dfrac{GroundTruth \cdot (1-GroundTruth)}
      #' {N_{Total}}
      #' }
      #'
      #' Where:
      #' \enumerate{
      #'   \item \eqn{R}: casual lucky guess.
      #'   \item \eqn{ProdAcc}: producer accuracy.
      #'   }
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #' @return A list with Ground Truth indexes, their variance, confidence
      #' intervals and the matrix with the expected frequencies.
      #'
      #' @examples
      #' A<-matrix(c(148,1,8,2,0,0,50,15,3,0,1,6,39,7,1,1,0,6,25,1,1,0,0,1,6),nrow=5,
      #' ncol=5)
      #' p<-ConfMatrix$new(A,Source="Türk 1979")
      #' p$GroundTruth()
      #'
      #' @aliases NULL


    GroundTruth=function(a=NULL){
      M<-self$Values
      M_0<-M-diag(diag(M),nrow(M),nrow(M))
      fi<-apply(M_0,1,sum)/sum(M_0)
      fj<-apply(M_0,2,sum)/sum(M_0)

      V<-list()
      U<-list()
      vi<-c()
      uj<-c()
      vi<-fi
      uj<-fj
      V[[1]]<-vi
      U[[1]]<-uj
      tol<-1
      i=1
      while (tol>=0.0001) {
        if((i)%%2==0){
          U[[i+1]]<-U[[1]]/(sum(V[[i]])-V[[i]])
        }else{
          V[[i+1]]<-V[[1]]/(sum(U[[i]])-U[[i]])
        }
        if(i>=3){
          if(!is.null(V[[i-1]])){
            tol<-abs(sum(V[[i+1]])-sum(V[[i-1]]))
          }else{
            tol<-abs(sum(U[[i+1]])-sum(U[[i-1]]))
          }
        }
        i=i+1
      }
      fij<-U[[length(U)]]%*%t(V[[length(V)]])
      fij_0<-fij-diag(diag(fij))
      Expfij<-fij_0*sum(M_0)
      Ri<-c()
      for (i in 1:length(V[[length(V)]])) {
        Ri<-c(Ri,V[[length(V)]][i]/sum(V[[length(V)]]))

      }
      Ai<-self$ProdAcc()[[1]]
      GroundTruth<-(Ai-Ri)/(1-Ri)
      VarGroundTruth<-(GroundTruth*(1-GroundTruth))/sum(self$Values)
      ConfInt<-list()
      for (i in 1:length(GroundTruth)) {
        ConfInt[[i]]<-c(private$ConfInt(GroundTruth[i],VarGroundTruth[i],
        a)$ConfInt_inf,private$ConfInt(GroundTruth[i],VarGroundTruth[i],
        a)$ConfInt_sup)
      }

    return(list(GroundTruth=GroundTruth,VarGroundTruth=VarGroundTruth,
                Conf_Int=ConfInt,ExpFrec=Expfij))
    },



      #' @description Public method that calculates the Ground Truth index
      #' for class i, its variance and confidence interval.The reference
      #' \insertCite{turk1979gt;textual}{ConfMatrix} is followed for the computations.
      #'
      #' To calculate R_i we begin the following iterative process:
      #' Be \eqn{U_j^{(0)}=f_j^0} with \eqn{f_j^0=\dfrac{\overline{x}_{i+}}{\sum_{i=1}^M \overline{x}_{i+}}}
      #' and \eqn{f_i^0=\dfrac{\overline{x}_{+i}}{\sum_{i=1}^M \overline{x}_{+i}}}
      #'
      #'  Where \eqn{2m} with \eqn{m=1,2,\cdots}
      #'  \deqn{
      #'  V_{i,2m-1}=\dfrac{f_i^0}{U_{+,2m-2}-U_{i,2m-2}}
      #'  }
      #'  where \eqn{U_{+,2m}=\sum_{i=1}^M U_{j,2m} }
      #'  and when \eqn{2m+1} with \eqn{m=1,2,\cdots}
      #'  \deqn{
      #'  U_{j,2m}=\dfrac{f_j^0}{V_{+,2m-1}-V_{i,2m-1}}
      #'  }
      #'  where \eqn{V_{+,2m-1}=\sum_{i=1}^M V_{i,2m-1} }
      #'
      #'  The iterative steps continue for \eqn{m=1, 2,\cdots} until
      #'  the accuracy stabilizes thus taking the V term.
      #'  Where
      #'  \deqn{
      #' R_i=\dfrac{V_i}{\sum_{i=1}^{k} V_i}
      #' }
      #'
      #' \deqn{
      #' ProdAcc_i=\dfrac{x_{ii}}{\sum_{j=1}^M x_{+j}}
      #' }
      #' \deqn{
      #' GroundTruth_i = \dfrac{ProdAcc_i-R_i}{1-R_i}
      #' }
      #'
      #' \deqn{
      #' \sigma^2_{GroundTruth_i}=\dfrac{GroundTruth_i \cdot (1-GroundTruth_i)}
      #' {N_{Total}}
      #' }
      #'
      #' Where:
      #' \enumerate{
      #'   \item \eqn{R_i}: casual lucky guess for class \eqn{i}. Is a real value.
      #'   \item \eqn{ProdAcc_i}: producer accuracy for class \eqn{i}.
      #'   }
      #' @param i \verb{
      #' Class to evaluate, where} \eqn{i \in \mathbb{Z}-\{0\}}.
      #'
      #' @param a \verb{
      #' Significance level. By default 0.05.
      #' }
      #' @return A list with Ground Truth index for class \eqn{i}, its variance, confidence
      #' interval and the matrix with the expected frequencies for all classes.
      #'
      #' @examples
      #' A<-matrix(c(148,1,8,2,0,0,50,15,3,0,1,6,39,7,1,1,0,6,25,1,1,0,0,1,6),nrow=5,
      #' ncol=5)
      #' p<-ConfMatrix$new(A,Source="Türk 1979")
      #' p$GroundTruth_i(3)
      #'
      #' @aliases NULL


    GroundTruth_i=function(i,a=NULL){
     GroundTruth<-self$GroundTruth()[[1]][i]
     VarGroundTruth<-self$GroundTruth()[[2]][i]
     ConfInt<-self$GroundTruth()[[3]][[i]]
     Expfij<-self$GroundTruth()[[4]]
    return(list(GroundTruth=GroundTruth,VarGroundTruth=VarGroundTruth,
                Conf_Int=ConfInt,ExpFrec=Expfij))
    },



      #' @description Public method that provides that Hellinger distance 
      #' between two confusion matrices.
      #' The reference \insertCite{garcia2018;textual}{ConfMatrix} is followed
      #' for the computations.
      #'
      #' \deqn{
      #' HellingerDist = \dfrac{4n_{A}m_{B}}{n_{A}+m_{B}} \sum^{M}_{i=1} (\sqrt{p_i}-\sqrt{q_i})^2
      #' }
      #'
      #' Where:
      #' \enumerate{
      #'   \item \eqn{n_{A}}: sum of elements of the matrix A.
      #'   \item \eqn{m_{B}}: sum of elements of the matrix B.
      #'   \item \eqn{p_i}: probability that element \eqn{i \in [1, \cdots, MxM]} is well classified in matrix A.
      #'   \item \eqn{q_i}: probability that element \eqn{i \in [1, \cdots, MxM]} is well classified in matrix B.
      #' }
      #' @param f \verb{
      #' Element of the ConfMatrix.
      #' }
      #' @param p \verb{
      #' probability vector of matrix A. By default, relative frequencies observed
      #' for each cell is taken.
      #' }
      #' @param q \verb{
      #' probability vector of matrix B. By default, relative frequencies observed
      #' for each cell is taken.
      #' }
      #'
      #' @return A real value for the Hellinger distance.
      #' @examples
      #' A<-matrix(c(65,6,0,4,4,81,11,7,22,5,85,3,24,8,19,90),nrow=4,ncol=4)
      #' r<-ConfMatrix$new(A,Source="Congalton and Green 2008")
      #' B<-matrix(c(45,6,0,4,4,91,8,7,12,5,55,3,24,8,9,55),nrow=4,ncol=4)
      #' f<-ConfMatrix$new(B,Source="Congalton and Green 2008")
      #' r$HellingerDist(f)
      #'
      #' @aliases NULL

    HellingerDist = function(f,p=NULL,q=NULL){
      if(class(f)[1]!="ConfMatrix"){
       stop("A ConfMatrix element is not being introduced\n")
      }
      A<-self$Values
      B<-f$Values
      if(is.null(p)){
      p<-A/sum(A)
      }else{p<-p}
      if(is.null(q)){
        q<-B/sum(B)
      }else{q<-q}

      if(length(q)!=length(p)){
        stop("Probabilities with different sizes.")
      }else{
        p_orig <- 4*((sum(self$Values)*sum(B)/(sum(self$Values)+sum(B))) *
                       sum((sqrt(p) - sqrt(q))^2))

      return(HellingerDist=p_orig)
      }
    },


# Functions that return multiple indices ----------------------------------


      #' @description Public method that calculates the values of quantity
      #' difference, exchange and shift. Quantity difference is the amount of
      #' difference between the product and the reference and is due to the
      #' less than maximum match in the proportions of the categories. Exchange
      #' represents transitions from class \eqn{i} to \eqn{j} and a transition from class \eqn{j}
      #' to class \eqn{i} in an identical number of cases. Shift refers to the
      #' difference remaining after subtracting quantity difference and exchange
      #' from the overall difference. The reference
      #' \insertCite{pontius2014;textual}{ConfMatrix} is followed for the computations.
      #'
      #' Where
      #' \deqn{
      #' Q=\dfrac{\sum^M_{j=1} q_{j}}{2}
      #'}
      #'
      #' \deqn{
      #' E=\dfrac{\sum^M_{j=1} e_{j}}{2}
      #'}
      #'
      #' \deqn{
      #' S=\dfrac{\sum^M_{j=1} s_{j}}{2}
      #' }
      #'
      #' with
      #'
      #' \deqn{
      #' d_{j}=\dfrac{ \left( \sum^M_{i=1} (x_{ij} + x_{ji})  \right) -2 x_{jj} }{\sum^M_{i=1} \sum^M_{j=1} x_{ij}}
      #'}
      #'
      #' \deqn{
      #' q_{j}=\dfrac{\left|  \sum^M_{i=1} (x_{ij} + x_{ji})  \right| }{\sum^J_{i=1} \sum^J_{j=1} x_{ij}}
      #'}
      #'
      #' \deqn{
      #' e_{j}=\dfrac{2 \left( \left( \sum^M_{i=1} min(x_{ij}, x_{ji})  \right) - x_{jj} \right)}{\sum^M_{i=1} \sum^M_{j=1} x_{ij}}
      #'}
      #'
      #' \deqn{
      #' s_{j}=d_{j}-q_{j}-e_{j}
      #'}
      #'
      #' @return A list of integer values with quantity, exchange, and shift.
      #' In addition to the differences for classes of the components of
      #' quantity, exchange and turn.
      #' @examples
      #' A<-matrix(c(3,2,1,1,3,3,2,0,1),nrow=3,ncol=3)
      #' p<-ConfMatrix$new(A,Source="Pontius Jr. and Santacruz 2023")
      #' p$QES()
      #'
      #' @aliases NULL

     QES = function(){
      nc <- nrow(self$Values)
      SumaMatriz <-sum(self$Values)
      SumaDigonal<-sum(diag(self$Values))
      ee<- matrix(rep(0, nc*nc), nrow =nc, ncol=nc, byrow=TRUE)
      d<- matrix(rep(0, nc), nrow =1, ncol=nc, byrow=TRUE)
      q<- matrix(rep(0, nc), nrow =1, ncol=nc, byrow=TRUE)
      e<- matrix(rep(0, nc), nrow =1, ncol=nc, byrow=TRUE)
      s<- matrix(rep(0, nc), nrow =1, ncol=nc, byrow=TRUE)
        for (j in 1:nc){
          for (i in 1:nc){
            if (i>j){
            ee[i,j] <-  (min(self$Values[i,j], self$Values[j,i]))*2
            }else{
              ee[i,j] <- 0
            }
          }
        }

        for (j in 1:nc){
        d[j]<-d[j]+ sum(self$Values[,j])+sum(self$Values[j,])-2*self$Values[j,j]
        q[j]<-q[j]+ abs(sum(self$Values[,j])- sum(self$Values[j,]))
        e[j]<-e[j]+ sum(ee[,j])+sum(ee[j,])
        s[j]<-d[j]- q[j]-e[j]
        }

      D<- sum(d)/2
      Q<- sum(q)/2
      E<- sum(e)/2
      S<- sum(s)/2

     return(list( OverallQuant=Q, OverallExch=E, OverallShift=S,
                  quant=q, exch=e, shift=s))
     },



# Functions that return matrices ------------------------------------------

      #' @description Public method that typifies the confusion matrix.
      #' The total sum of the original matrix is used for typing. In a
      #' typed matrix the sum of all values is unity. The resulting
      #' values can be presented as real values (parameter RaR=1), or as
      #' a percentage (parameter RaR !=1).
      #'
      #'  \deqn{
      #' MTypify=\dfrac{x_{ij}}{\sum^M_{i,j=1} x_{ij}}
      #' }
      #'
      #' @param RaR "1" indicates result as real, other values mean percentage
      #' as integer. By default RaR=1.
      #' @return A list with two arrays, the first is the original array,
      #' the second the typed one.
      #' @examples
      #' A<-matrix(c(65,6,0,4,4,81,11,7,22,5,85,3,24,8,19,90),nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A, Source="Congalton and Green 2008")
      #' p$MTypify(RaR=5)

     MTypify =function(RaR=NULL){
        if(!is.null(RaR)){
          RaR <- RaR
        }else{R<-1}
       
      ConfM=self$Values
      M <- ConfM/(sum(ConfM))
        if (RaR==1){
        M <- ConfM/(sum(ConfM))
        return(M)
        } else {
      M <- ConfM/(sum(ConfM))
      M[] <- as.integer(100*M)
      return(list(OriginalMatrix=self$Values,TypifyMatrix=M))
      }
     },



      #' @description Public method that provides B resamples, using a
      #' multinomial distribution, of the confusion matrix of a ConfMatrix
      #' object. As a result, a set of bootstrapped cases is offered. The
      #' reference \insertCite{fienberg1970;textual}{ConfMatrix} is
      #' followed for the computations.
      #' @param B Number of resamples.
      #' @param pr Vector with resampling probabilities. By default, the
      #' success probability of each cell will be taken.
      #' @return A list of B + 1 arrays formed by the original confusion matrix
      #' and all the simulated cases.
      #' @examples
      #' A<-matrix(c(65,6,0,4,4,81,11,7,22,5,85,3,24,8,19,90),nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A, Source="Congalton and Green 2008")
      #' p$MBootStrap(2)
      #'
      #' @aliases NULL

     MBootStrap=function(B,pr=NULL){
      nc<-ncol(self$Values)
      M1<-self$Values
      if(is.null(pr)){
        pr<-M1/sum(M1)
      }else{
        pr<-pr
      }
      M2<-list()
      boots<-rmultinom(B,sum(M1),pr)
        for(i in 1:ncol(boots)){
          M2[[i]]<-matrix(boots[,i],ncol=nc,nrow=nc)
        }

     return(list(OriginalMatrix=self$Values,BootStrap=M2))
     },


      #' @description Public method that carries out an iterative process in
      #' order to equals one the sum of values by rows and columns.
      #' The references \insertCite{fienberg1970;textual}{ConfMatrix} and
      #'  \insertCite{munoz2016;textual}{ConfMatrix} are followed for the computations.
      #'
      #' The following iterative process is used:
      #'
      #' Let \eqn{x_{ij}} be the elements of the instance. It defines:
      #'
      #' \deqn{x'_{ij}=\dfrac{x_{ij}}{x_{i+}}}
      #'
      #'\deqn{x''_{ij}=\dfrac{x'_{ij}}{x'_{+j}}}
      #'
      #' Taking \eqn{x_{ij}=x''_{ij}} for the next iteration.
      #'
      #' @param iter \verb{
      #' Number of iteration. By default iter=1000.
      #' }
      #' @return A list formed by the original confusion matrix and the
      #' normalized matrix.
      #'
      #' @examples
      #' A<-matrix(c(238051,7,132,0,0,24,9,2,189,1,4086,188,0,4,16,45,1,0,939,5082,
      #' 51817,0,34,500,1867,325,17,0,0,5,11148,1618,78,0,0,0,0,48,4,834,2853,340,
      #' 32,0,197,5,151,119,135,726,6774,75,1,553,0,105,601,110,174,155,8257,8,0,
      #' 29,36,280,0,0,6,5,2993,0,115,2,0,4,124,595,0,0,4374),nrow=9,ncol=9)
      #' p<-ConfMatrix$new(A,Source="Muñoz 2016")
      #' p$MNormalize()
      #'
      #' @aliases NULL

     MNormalize=function(iter=NULL){

      if(!is.null(iter)){
        iter <- iter
      }else{iter<-1000}

      rg<-nrow(self$Values)
      x1<-self$Values
        for (k in 1:iter) {
        sumfilas=apply(x1,1,sum)
          for (i in 1:rg) {
          x1[i,]=x1[i,]/sumfilas[i]
          }
        sumcolumnas=apply(x1,2,sum)
          for (j in 1:rg) {
          x1[,j]=x1[,j]/sumcolumnas[j]
          }
        }
      NormMatrix<-x1
     return(list(OriginalMatrix=self$Values,NormalizeMatrix=NormMatrix))
     },



      #' @description Public method that small values are calculated for empty
      #' cells of the matrix. All non-empty cells of the matrix change their
      #' values. This function will not be applied if all the elements of the
      #' matrix are different from 0.
      #' The reference \insertCite{munoz2016;textual}{ConfMatrix} is followed
      #' for the computations.
      #'
      #' Let \eqn{x_{ij}} be the elements of the instance.
      #'
      #' The following values are defined:
      #'
      #' \deqn{e_{ij}=\dfrac{x_{i+}x_{+j}}{\sum^M_{i,j=1} x_{ij}}}
      #'
      #' \deqn{v=\dfrac{\left( \sum^M_{i,j=1} x_{ij} \right)^2 - \sum^M_{i,j=1} x_{ij}^2}{\sum^M_{i,j=1} (e_{ij}-x_{ij})^2}}
      #'
      #' \deqn{p_{ij}=\dfrac{e_{ij} \cdot v }{\sum^M_{i,j=1} x_{ij}}}
      #'
      #' Finally, the elements of the pseudozero matrix \eqn{Z} will be given by:
      #'
      #' \deqn{z_{ij}=\left(\dfrac{\sum^M_{i,j=1} x_{ij}}{(\sum^M_{i,j=1} x_{ij})+v} \right)
      #' (p_{ij}+x_{ij})}
      #'
      #' @return A list formed by the original confusion matrix and the
      #' Pseudozeroes matrix.
      #' @examples
      #' A<-matrix(c(238051,7,132,0,0,24,9,2,189,1,4086,188,0,4,16,45,1,0,939,5082,
      #' 51817,0,34,500,1867,325,17,0,0,5,11148,1618,78,0,0,0,0,48,4,834,2853,340,
      #' 32,0,197,5,151,119,135,726,6774,75,1,553,0,105,601,110,174,155,8257,8,0,
      #' 29,36,280,0,0,6,5,2993,0,115,2,0,4,124,595,0,0,4374),nrow=9,ncol=9)
      #' p<-ConfMatrix$new(A,Source="Muñoz 2016")
      #' p$MPseudoZeroes()
      #'
      #' @aliases NULL

     MPseudoZeroes = function(){
       k=0
       rg<-nrow(self$Values)
        for (i in 1:rg) {
         for (j in 1:rg) {
          if((self$Values[i,j]!=0)==TRUE){
            k=k+1
             if(k==length(self$Values)){
               stop("\nThe Pseudoceros Matrix removes the zeros from the matrix.
                    \nYour matrix does not have any zeros to remove.\n")
             }}
        }
       }

       ConfM=self$Values
       SumaMatriz <-sum(ConfM)
       MLandas <- (private$sumfil(self$Values) %*% t(private$sumcol(self$Values)))/(SumaMatriz*SumaMatriz)
       K <- (SumaMatriz*SumaMatriz -
               sum(ConfM*ConfM))/sum((SumaMatriz*MLandas - ConfM )^2)
       MPseudoceros <- (SumaMatriz/(K+SumaMatriz))*(ConfM + K*MLandas)

     return(list(OriginalMatrix=self$Values,PseudoZeroesMatrix=MPseudoceros))
     },


# test function -----------------------------------------------------------



      #' @description Public method that tests whether two independent
      #' confusion matrices (instances of the ConfMatrix class), are
      #' significantly different using their overall accuracy indexes.
      #' The reference \insertCite{congalton2008;textual}{ConfMatrix} and \insertCite{ma1995Tau;textual}{ConfMatrix} are followed
      #' for the computations.
      #'
      #' \deqn{
      #' Z = \dfrac{|O_A-O_B|}{\sqrt{(\sigma^2_{O_A}+\sigma^2_{O_B})}}
      #' }
      #'
      #' Where:
      #' \enumerate{
      #'   \item \eqn{O_A}: overall index of matrix A.
      #'   \item \eqn{O_B}: overall index of matrix B.
      #'   \item \eqn{\sigma^2_{O_A}}: variance of \eqn{O_A}.
      #'   \item \eqn{\sigma^2_{O_B}}: variance of \eqn{O_B}.
      #' }
      #' @param f \verb{
      #' Instance of ConfMatrix class.
      #' }
      #' @return A list of class "htest" containing the results of the hypothesis test.
      #'
      #' @examples
      #' A<-matrix(c(65,6,0,4,4,81,11,7,22,5,85,3,24,8,19,90),nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A,Source="Congalton and Green 2008")
      #' B<-matrix(c(45,6,0,4,4,91,8,7,12,5,55,3,24,8,9,55),nrow=4,ncol=4)
      #' f<-ConfMatrix$new(B,Source="Congalton and Green 2008")
      #' p$OverallAcc.test(f)
      #'
      #' @aliases NULL

    OverallAcc.test=function(f){
      if(class(f)[1]!="ConfMatrix"){
        stop("A ConfMatrix element is not being introduced")
      }

      k1<-self$OverallAcc()[[1]]
      k2<-f$OverallAcc()[[1]]
      v1<-self$OverallAcc()[[2]]
      v2<-f$OverallAcc()[[2]]

      Z<-abs(k1-k2)/sqrt(v1+v2)

      htest_result <- list(
        statistic = c(Z=round(Z,4)),    
        p.value = 2 * (1 - pnorm(abs(Z))),
        method = "Test for overall accuracy difference",
        data.name = paste("Overall accuracy of", self$ID ,"vs overall accuracy of",f$ID) 
      )

      class(htest_result) <- "htest"
      return(htest_result)
    },


      #' @description Public method that tests whether two independent
      #' confusion matrices (instances of the ConfMatrix class), are
      #' significantly different when using the kappa indexes.
      #' The reference \insertCite{congalton2008;textual}{ConfMatrix} is followed
      #' for the computations.
      #'
      #' \deqn{
      #' Z = \dfrac{|k_A-k_B|}{\sqrt{(\sigma^2_{k_A}+\sigma^2_{k_B})}}
      #' }
      #'
      #' Where:
      #' \enumerate{
      #'   \item \eqn{k_A}: kappa index of matrix A.
      #'   \item \eqn{k_B}: kappa index of matrix B.
      #'   \item \eqn{\sigma^2_{k_A}}: variance of \eqn{k_A}.
      #'   \item \eqn{\sigma^2_{k_B}}: variance of \eqn{k_B}.
      #' }
      #' @param f \verb{
      #' Element of the ConfMatrix class.
      #' }
      #' @return A list of class "htest" containing the results of the hypothesis test.
      #'
      #' @examples
      #' A<-matrix(c(65,6,0,4,4,81,11,7,22,5,85,3,24,8,19,90),nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A,Source="Congalton and Green 2008")
      #' B<-matrix(c(45,6,0,4,4,91,8,7,12,5,55,3,24,8,9,55),nrow=4,ncol=4)
      #' f<-ConfMatrix$new(B,Source="Congalton and Green 2008")
      #' p$Kappa.test(f)
      #'
      #' @aliases NULL

      Kappa.test=function(f){
        if(class(f)[1]!="ConfMatrix"){
          stop("A ConfMatrix element is not being introduced")
        }

        k1<-self$Kappa()[[1]]
        k2<-f$Kappa()[[1]]
        v1<-self$Kappa()[[2]]
        v2<-f$Kappa()[[2]]

        Z<-abs(k1-k2)/(sqrt(v1+v2))

        htest_result <- list(
          statistic = c(Z=round(Z,4)),        
          p.value = 2 * (1 - pnorm(abs(Z))),
          method = "Test for kappa difference", 
          data.name = paste("Kappa index of", self$ID ,"vs Kappa index of", f$ID) 
        )

        class(htest_result) <- "htest"
        return(htest_result)
      },



      #' @description Public method that tests whether two independent
      #' confusion matrices (instances of the ConfMatrix class), are
      #' significantly different using their Tau indexes.
      #' The reference \insertCite{congalton2008;textual}{ConfMatrix} and
      #' \insertCite{ma1995Tau;textual}{ConfMatrix} are followed for the computations.
      #'
      #' \deqn{
      #' Z = \dfrac{|\tau_A-\tau_B|}{\sqrt{(\sigma^2_{\tau_A}+\sigma^2_{\tau_B})}}
      #' }
      #'
      #' Where:
      #' \enumerate{
      #'   \item \eqn{\tau_A}: Tau index of matrix A.
      #'   \item \eqn{\tau_B}: Tau index of matrix B.
      #'   \item \eqn{\sigma^2_{\tau_A}}: variance of \eqn{\tau_A}.
      #'   \item \eqn{\sigma^2_{\tau_B}}: variance of \eqn{\tau_B}.
      #' }
      #' @param f \verb{
      #' Element of the ConfMatrix class.
      #' }
      #' @return A list of class "htest" containing the results of the hypothesis test.
      #'
      #' @examples
      #' A<-matrix(c(65,6,0,4,4,81,11,7,22,5,85,3,24,8,19,90),nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A,Source="Congalton and Green 2008")
      #' B<-matrix(c(45,6,0,4,4,91,8,7,12,5,55,3,24,8,9,55),nrow=4,ncol=4)
      #' f<-ConfMatrix$new(B,Source="Congalton and Green 2008")
      #' p$Tau.test(f)
      #'
      #' @aliases NULL

    Tau.test=function(f){
      if(class(f)[1]!="ConfMatrix"){
        stop("A ConfMatrix element is not being introduced")
      }

      k1<-self$Tau()[[1]]
      k2<-f$Tau()[[1]]
      v1<-self$Tau()[[2]]
      v2<-f$Tau()[[2]]

      Z<-abs(k1-k2)/sqrt(v1+v2)

      htest_result <- list(
        statistic = c(Z=round(Z,4)),                    
        p.value = 2 * (1 - pnorm(abs(Z))),
        method = "Test for tau difference", 
        data.name = paste("Tau index of",self$ID,"vs Tau index of",f$ID) 
      )

      class(htest_result) <- "htest"
      return(htest_result)
    },



      #' @description Public method that performs a homogeneity test
      #' based on the Hellinger distance between two confusion matrices
      #' (instances of the ConfMatrix class).
      #' The test considers the individual cell values in the matrices.
      #' Bootstrap is applied to the matrices to obtain a consistent estimator.
      #' The reference \insertCite{garcia2018;textual}{ConfMatrix} are followed for
      #' the computations.
      #' The calculation consists of obtaining a statistic, which we will call
      #' \eqn{T_{n,m}}, between both matrices from ConfMatrix$HellingerDist.
      #' Bootstrap is then applied to the confusion matrices to obtain
      #' simulations of both matrices. ConfMatrix$HellingerDist is applied
      #' again to these simulations and we will obtain the statistics
      #' \eqn{T^*_{n,m}}. The p value is defined as:
      #' \deqn{
      #' \hat{p}=\dfrac{Card(T^*_{n,m} \geq T^*_{n,m})}{B}
      #' }
      #' @param B \verb{
      #' Number of bootstraps that you want to generate. By default B=1000.
      #' }
      #' @param f \verb{
      #' Element of the ConfMatrix class.
      #' }
      #' @return A list of class "htest" containing the results of the hypothesis test.
      #'
      #' @examples
      #' A<-matrix(c(65,6,0,4,4,81,11,7,22,5,85,3,24,8,19,90),nrow=4,ncol=4)
      #' p<-ConfMatrix$new(A,Source="Congalton and Green 2008")
      #' C<-matrix(c(45,6,0,4,4,91,8,7,12,5,55,3,24,8,9,55),nrow=4,ncol=4)
      #' f<-ConfMatrix$new(C,Source="Congalton and Green 2008")
      #' p$TSCM.test(f)
      #'
      #' @aliases NULL

    TSCM.test=function(f,B=NULL){
      if(class(f)[1]!="ConfMatrix"){
        stop("\nA ConfMatrix element is not being introduced.\n")
      }
      if(is.null(B)){
        B<-1000
      }else{B<-B}

      A<-self$Values
      C<-f$Values
      n<-length(A)
      m<-length(C)
      p2<-A/sum(A)
      q2<-C/sum(C)
      p_orig<-self$HellingerDist(f)
      p_0<-c()
      for(i in 1:length(p2)){
        p_01<-(n*p2[i]+m*q2[i])/(n+m)
        p_0<-c(p_0,p_01)
      }

      q1<-list()
      p1<-list()
      p1<-self$MBootStrap(B,p_0)[[2]]
      q1<-f$MBootStrap(B,p_0)[[2]]
      Tn<-c()

      for (i in 1:length(p1)) {
        Tn1<-self$HellingerDist(f,p=(p1[[i]]/sum(p1[[i]])),q=(q1[[i]]/sum(q1[[i]])))
        Tn<-c(Tn,Tn1)
      }

      Tn_boot<-c()
      for (i in 1:length(Tn)) {
        if((Tn[i]>=p_orig)==TRUE){
          Tn_boot<-c(Tn_boot,Tn[i])
        }
      }
      pvalue<-length(Tn_boot)/length(Tn)

      htest_result <- list(
        p.value = pvalue,        
        method = "TSCM-test for evaluate the similarity between two confusion matrices.\nBased on the discrete squared
        Hellinger distance", 
        data.name = paste(self$ID , "vs", f$ID) 
      )

      class(htest_result) <- "htest"
      return(htest_result)
    },

      #' @description Public method that performs the
      #'  quasi-independence test for the elements of a confusion matrix.
      #'  The reference
      #'  \insertCite{turk1979gt;textual}{ConfMatrix} and \insertCite{goodman1968analysis;textual}{ConfMatrix}
      #'  are followed for the computations.
      #'
      #' \deqn{
      #' G^2 = 2 \cdot \sum \log \dfrac{x_{ij}}{E_{ij}}
      #' }
      #'
      #' Following the procedure for calculating the elements
      #' of the function ConfMatrix$GroundTruth, we will have to \eqn{E_{ij}}
      #' is obtained from:
      #'
      #' \deqn{f_{ij}=U_j \cdot V_i}
      #' \deqn{f_{ij}^0=f_{ij}-f_{ii}}
      #' \deqn{M^0=x_{ij}-x_{ii}}
      #' where the elements of \eqn{M^0} are \eqn{m_{ij}^0}
      #'
      #' \deqn{E_{ij}=f_{ij}^0 \sum^M_{i,j=1} m_{ij}^{0}}
      #'
      #'
      #'
      #' Where:
      #' \enumerate{
      #'   \item \eqn{x_{ij}}: matrix element. Observed frequency.
      #'   \item \eqn{E_{ij}}: expected frequency.
      #' }
      #' @return A list of class "htest" containing the results of the hypothesis test.
      #'
      #' @examples
      #' A<-matrix(c(148,1,8,2,0,0,50,15,3,0,1,6,39,7,1,1,0,6,25,1,1,0,0,1,6),nrow=5,
      #' ncol=5)
      #' p<-ConfMatrix$new(A,Source= "Türk 1979")
      #' p$QIndep.test()
      #'
      #' @aliases NULL


    QIndep.test=function(){
       A_0<-self$Values-diag(diag(self$Values),nrow(self$Values),
          nrow(self$Values))
      Expfij<-self$GroundTruth()[[4]]
      k<-ncol(A_0)*ncol(A_0)-3*ncol(A_0)+1
      matr2<-A_0/Expfij
      matr2[is.nan(matr2)] <- 0

      for (i in 1:nrow(A_0)) {
        for (j in 1:nrow(A_0)) {
          if(matr2[i,j]==0){
            matr2[i,j]<-0
          }else{
            matr2[i,j]<-log(matr2[i,j])
          }
        }
      }

      Z<-2*sum(matr2)
      p.value <- pchisq(Z, df = k, lower.tail = FALSE)

      htest_result <- list(
        statistic = c(Z=round(Z,4)),     
        p.value = p.value,
        method = "Quasi-Independence Test",
        data.name = paste("ConfMatrix")
      )

      class(htest_result) <- "htest"
      return(htest_result)

    }

   ),

# Private functions -------------------------------------------------------

   private = list(

     sequence = function() {
       variables <- ls(envir = .GlobalEnv)
       es_confmatrix <- sapply(variables, function(x) {
         obj <- tryCatch(get(x, envir = .GlobalEnv), error = function(e) NULL)
         if (!is.null(obj)) {
           return(inherits(obj, "ConfMatrix"))
         } else {
           return(FALSE)
         }
       }, USE.NAMES = FALSE)
       
       number_confmatrix <- sum(es_confmatrix, na.rm = TRUE)
       return(number_confmatrix + 1)
     },

     
     ConfInt=function(p,var,a=NULL){
       if(is.null(a)){
         a<-0.05
       }else{a<-a}
       z<-qnorm(1-a/2)
       ConfInt_inf<-p-z*sqrt(var)
       ConfInt_sup<-p+z*sqrt(var)
       return(list(ConfInt_inf=ConfInt_inf,ConfInt_sup=ConfInt_sup))
     },

     sumfil=function(v){
       sumfil<-apply(v,1,sum)
       return(sumfil)
     },

     sumcol=function(v){
       sumcol<-apply(v,2,sum)
       return(sumcol)
     }

   ),
   active = list(
   )
)
