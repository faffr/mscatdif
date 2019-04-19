#' Multistage CAT DIF Function
#'
#' This function performs DIF analysis on items from a multi-stage CAT
#' @param scored.resp is a matrix of scored responses.
#' @param person.ability is a vector of person ability estimates.
#' @param item.difficulty is a vector of item difficulty measures.
#' @param group.info is a vector of person demographic variable.
#' @param which.focal is the desired focal group. Random group chosen if NULL.
#' @param which.ref is the desired reference group. Random group chosen if NULL.
#' @param item.fit is a vector of item fit information. Default = NULL.
#' @param item.names is a vector of item names. Default = NULL. 
#' @param outfilename is the name of the desired .csv output file. Default = NULL.
#' @param purification run purification? Default = FALSE.
#' @param field.test is it a field test? Default = FALSE.
#' @param FTnSize if !is.null(field.test) then FTnSize integer scalar must be specified.  
#' @keywords DIF multistage CAT
#' @export
#' @examples
#' set.seed(1010)
#' ItemData = matrix(sample(rep(0:1, each = 1000),1000),100,10)
#' PersMeas = rnorm(100, 0, 1)
#' DemogGroup = sample(rep(c("O","H"), each = 100), 100)
#' ItemMeas = rnorm(10, 0, 1)
#' mscatdif(ItemData, PersMeas, ItemMeas, DemogGroup)

mscatdif = function(scored.resp, person.ability, item.difficulty, group.info, which.focal = NULL, which.ref = NULL, item.fit = NULL, item.names = NULL, outfilename = NULL, purification = FALSE, field.test = FALSE, FTnSize = NULL){

  if(is.null(which.focal)){which.focal = names(table(group.info) )[1]} #whichever comes first in alphabet
  if(is.null(which.ref)){which.ref = names(table(group.info) )[2]}
  if(!is.null(outfilename)){outfilename = paste0(outfilename,".csv")}
  if(is.null(item.names)){item.names = paste0("Item ", seq(1:length(item.difficulty))) }
  if(field.test == TRUE){FTdifficulty = item.difficulty[1:FTnSize]}
  if(field.test == TRUE){item.difficulty = item.difficulty[(FTnSize+1):length(item.difficulty)]}

  #CREATE A DATA.FRAME FOR ALL RESP, ABILITY, AND GROUP
  alldata <- data.frame(scored.resp, measure = person.ability, group.info = group.info )

  alldata[,"group.info"] <- as.character(alldata[,"group.info"])

  #RENAME THE SCORED RESPONSE COLUMNS
  colnames(alldata)[1:ncol(scored.resp)] = paste0( "X" , seq(1,ncol(scored.resp)) )

  #CREATE THE EXPECTED SCORE VALUE FOR EVERY EXAMINEE
  alldata[,"expected"] = apply(sapply(item.difficulty,function(x) 1/(1+exp(-(person.ability-x)))),1,sum)

  #RELABEL GROUP BASED ON FOCAL.NAME
  #NOTE THAT ALL MISSING DATA MUST BE REMOVED PRIOR
  alldata$group.info[as.character(alldata$group.info) == which.focal] <- 1
  alldata$group.info[as.character(alldata$group.info) == which.ref] <- 0

  #NUMBER OF ITEMS
  its = length(item.difficulty)
  if(field.test == TRUE){its = FTnSize}

  #FINAL DATA FRAME WHERE MH STATISTICS IS TO BE COLLECTED

  MH_collection1 = data.frame(1:its)

  ######################################
  #APPLY STANDARD DIF METHOD
  ######################################

  #FUNCTION TO CARRY OUT MH DIF ANALYSIS

  MHfun <- function(MHdataframe){

    for(i in 1:its){

      #SUBSET INDIVIDUALS THAT HAVE RESPONSES FOR ITEM i
      nomissing = subset(alldata,!is.na(alldata[,i]))

      if(field.test == TRUE){
        #ADD FT ITEM i RESPONSE SCORE TO THE "expected" score
        nomissing[,"expected"] = nomissing[,"expected"] + nomissing[,i]
      }

      #CREATE GROUPS OF EXPECTED SCORES BY UNITS OF 2 FROM 0 TO #OF ITEMS
      nomissing[,"expected_cat"] = cut(nomissing[,"expected"], breaks = c(seq(0,max(nomissing[,"expected"]) ,1)) )


      #SUBSET FOCAL: FEMALES AND REFERENCE: MALES/ FOCAL: HISP AND REFERENCE: NON-HISPANICS AND CREATE A NUMERIC CATEGORY
      dataR = subset(nomissing,nomissing[,"group.info"]=="0")
      dataR[,"expected_cat_numeric"] = as.numeric(dataR[,"expected_cat"])

      dataF = subset(nomissing,nomissing[,"group.info"]=="1")
      dataF[,"expected_cat_numeric"] = as.numeric(dataF[,"expected_cat"])

      #GET THE MAX NUMBER OF CATEGORIES
      cats = max(dataR[,"expected_cat_numeric"],na.rm=TRUE)

      #CELL COUNTS ARE COLLECTED IN THIS MATRIX
      temp = matrix(rep(NA,8*cats),8,cats)

      #LOOP TO OBTAIN COUNTS FOR EACH SCORE CATEGORY
      for(j in 1:cats){


        #OBTAIN THE NUMBER OF INDIVIDUALS IN THE SAME GROUP AND THE SAME EXPECTED SCORE CATEGORY
        A=as.numeric( length(subset(dataR[,i],dataR[,"expected_cat_numeric"]==j & dataR[,i]==1)) )
        B=as.numeric( length(subset(dataR[,i],dataR[,"expected_cat_numeric"]==j & dataR[,i]==0)) )

        C=as.numeric( length(subset(dataF[,i],dataF[,"expected_cat_numeric"]==j & dataF[,i]==1)) )
        D=as.numeric( length(subset(dataF[,i],dataF[,"expected_cat_numeric"]==j & dataF[,i]==0)) )

        #GET THE DIFFERENT PARTS OF THE MH STATISTIC
        N = A+B+C+D
        temp[1,j] = A + D
        temp[2,j] = B + C
        temp[3,j] = A*D
        temp[4,j] = B*C
        temp[5,j] = N
        temp[6,j] = (A+B)*(A+C)/N
        temp[7,j] = A
        temp[8,j] = (A + B)*(C + D)*(A + C)*(B + D) / (N^2*(N-1))

      }

      #REMOVE ALL CELL COUNTS WITH ZERO
      temp2 = t(temp)
      temp2 = temp2[apply(temp2,1,function(y) !any(is.nan(y))),]
      temp = t(temp2)

      #CALCULATE THE ALPHA, X^2, MH-DIF, SEs, ETS CLASSIFICATION
      #ALPHA
      MHdataframe[i,"Alpha"] = sum(temp[3,]/temp[5,])/sum(temp[4,]/temp[5,])
      #CHI-SQUARE STATISTIC
      MHdataframe[i,"X^2"] = round(( ( abs(sum(temp[7,]) - sum(temp[6,])) - 0.5 )^2)/sum(temp[8,]),3)

      #P-VALUE
      MHdataframe[i,"P-Value"] = round(1-pchisq(MHdataframe[i,"X^2"],1),3)

      #MH D-DIFF
      MHdataframe[i,"MH D-DIF"] = -2.35*log(MHdataframe[i,"Alpha"])

      #SEs
      cats = ncol(temp)
      var = matrix(rep(NA,cats*2),2,cats)

      var[1,] = temp[3,] + MHdataframe[i,"Alpha"]*temp[4,]
      var[2,] = temp[1,] + MHdataframe[i,"Alpha"]*temp[2,]

      MHdataframe[i,"SE(MH D-DIF)"] = sqrt(sum(var[1,]*var[2,]/temp[5,]^2)/(2*(sum(temp[3,]/temp[5,]))^2))*2.35


      #ETS A,B,C RULE from Zwick 2012

      if( is.na(MHdataframe[i,"MH D-DIF"])  )

      { MHdataframe[i,"ETS Rule"] = NA  }

      else

        if( abs(MHdataframe[i,"MH D-DIF"]) >= 1 & (abs(MHdataframe[i,"MH D-DIF"])-1)/MHdataframe[i,"SE(MH D-DIF)"] > 1.645 &
            abs(MHdataframe[i,"MH D-DIF"]) >= 1.5 & MHdataframe[i,"MH D-DIF"] > 0 )

        { MHdataframe[i,"ETS Rule"] = "C"; MHdataframe[i,"Favored"] = which.focal }

      else

        if( abs(MHdataframe[i,"MH D-DIF"]) >= 1 & (abs(MHdataframe[i,"MH D-DIF"])-1)/MHdataframe[i,"SE(MH D-DIF)"] > 1.645 &
            abs(MHdataframe[i,"MH D-DIF"]) >= 1.5 & MHdataframe[i,"MH D-DIF"] < 0  )

        { MHdataframe[i,"ETS Rule"] = "C"; MHdataframe[i,"Favored"] = which.ref  }

      else

        if(  abs(MHdataframe[i,"MH D-DIF"]) >= 1 & abs(MHdataframe[i,"MH D-DIF"])/MHdataframe[i,"SE(MH D-DIF)"] > 1.96 &
             MHdataframe[i,"MH D-DIF"] > 0    )

        { MHdataframe[i,"ETS Rule"] = "B"; MHdataframe[i,"Favored"] = which.focal  }

      else

        if(  abs(MHdataframe[i,"MH D-DIF"]) >= 1 & abs(MHdataframe[i,"MH D-DIF"])/MHdataframe[i,"SE(MH D-DIF)"] > 1.96 &
             MHdataframe[i,"MH D-DIF"] < 0    )

        { MHdataframe[i,"ETS Rule"] = "B"; MHdataframe[i,"Favored"] = which.ref  }
      else

        if(  (abs(MHdataframe[i,"MH D-DIF"]) < 1 | MHdataframe[i,"P-Value"] > 0.05) & MHdataframe[i,"MH D-DIF"] > 0   )

        { MHdataframe[i,"ETS Rule"] = "A"; MHdataframe[i,"Favored"] = which.focal  }

      else

        if(  (abs(MHdataframe[i,"MH D-DIF"]) < 1 | MHdataframe[i,"P-Value"] > 0.05) & MHdataframe[i,"MH D-DIF"] < 0   )

        { MHdataframe[i,"ETS Rule"] = "A"; MHdataframe[i,"Favored"] = which.ref  }

      colnames(MHdataframe)[1] = "Item"

      #TOTAL ITEM LEVEL SAMPLE SIZE
      MHdataframe[i, paste0("N: ",which.ref)] = nrow(dataR)
      MHdataframe[i, paste0("N: ",which.focal)] = nrow(dataF)

      #ATTACH ITEM NAMES
      MHdataframe[i,"Names"] = item.names[MHdataframe$"Item"[i]]

      #PROPORTION CORRECT BY GROUP
      MHdataframe[i, paste0("Prop. ",which.ref)] = round(sum(as.numeric(dataR[,i] )/nrow(dataR) ), 3)
      MHdataframe[i, paste0("Prop. ",which.focal)] = round(sum(as.numeric(dataF[,i]) )/nrow(dataF) ,3)

      #ITEM DIFFICULTY
      if(field.test == TRUE){MHdataframe[i, "Difficulty"] = FTdifficulty[i]}else{MHdataframe[i, "Difficulty"] = item.difficulty[i]}


      #OVERALL PROP. CORRECT
      MHdataframe[i, "CTT pval"] = round( (sum(dataR[,i]) + sum(dataF[,i]) )/( nrow(dataR) + nrow(dataF) ), 3)

      MHdataframe[i,"Match Cat."] = cats

      print(paste0("item ",MHdataframe$"Item"[i]," - ",cats," expected score categories for matching"))


      #INFIT/OUTFIT STATISTICS
      if(!is.null(item.fit) ) {MHdataframe[i, "IN.MSQ"] = item.fit[i,1]}
      if(!is.null(item.fit) ) {MHdataframe[i, "OUT.MSQ"] = item.fit[i,2]}



    }

    return(MHdataframe)

  }

  #RUN THE FUNCTION

  MH_collection1 <- MHfun(MH_collection1)


  ######################################
  #APPLY EB DIF METHOD Zwick, 1997, 1999
  ######################################

  #FUNCTION TO RUN EB METHOD

  MHfunEB <- function(MHdataframe){

    mu = mean(MHdataframe[,"MH D-DIF"], na.rm = TRUE)

    tau2 = var(MHdataframe[,"MH D-DIF"], use = "complete.obs" ) - mean(MHdataframe[,"SE(MH D-DIF)"]^2, na.rm = TRUE )

    for(i in 1:its){

      MHdataframe[i,"W"] = tau2/(MHdataframe[i,"SE(MH D-DIF)"]^2+tau2)
      MHdataframe[i,"Post Mean"] = MHdataframe[i,"W"]*MHdataframe[i,"MH D-DIF"]+(1-MHdataframe[i,"W"])*mu
      MHdataframe[i,"Post SD"] = sqrt( (MHdataframe[i,"W"])*MHdataframe[i,"SE(MH D-DIF)"]^2  )
      MHdataframe[i,"Post pred SD"] = sqrt( (1 + MHdataframe[i,"W"])*MHdataframe[i,"SE(MH D-DIF)"]^2  )

      #TRUE DIF METHOD

      MHdataframe[i,"TRUE.C-"] = round( pnorm( -1.5,  MHdataframe[i,"Post Mean"], MHdataframe[i,"Post SD"] ) , 3)


      MHdataframe[i,"TRUE.B-"] = round( ( pnorm(-1.0, MHdataframe[i,"Post Mean"], MHdataframe[i,"Post SD"] ) - pnorm(-1.5, MHdataframe[i,"Post Mean"], MHdataframe[i,"Post SD"] ) )  , 3)


      MHdataframe[i,"TRUE.A"]  = round( pnorm(1, mean = MHdataframe[i,"Post Mean"], sd =  MHdataframe[i,"Post SD"] ) - pnorm(-1, mean = MHdataframe[i,"Post Mean"], sd =  MHdataframe[i,"Post SD"] ) , 3)


      MHdataframe[i,"TRUE.B+"] = round( ( pnorm( 1.5,  MHdataframe[i,"Post Mean"], MHdataframe[i,"Post SD"] ) - pnorm(1.0, MHdataframe[i,"Post Mean"], MHdataframe[i,"Post SD"] ) )  , 3)


      MHdataframe[i,"TRUE.C+"] = round( ( pnorm( 100,  MHdataframe[i,"Post Mean"], MHdataframe[i,"Post SD"] ) - pnorm(1.5, MHdataframe[i,"Post Mean"], MHdataframe[i,"Post SD"] ) ) , 3)



      #FUTURE DIF METHOD


      MHdataframe[i,"FUT.C-"] = round( pnorm( -1.5,  MHdataframe[i,"Post Mean"], MHdataframe[i,"Post pred SD"] ) , 3)


      MHdataframe[i,"FUT.B-"] = round( ( pnorm(-1.0, MHdataframe[i,"Post Mean"], MHdataframe[i,"Post pred SD"] ) - pnorm(-1.5, MHdataframe[i,"Post Mean"], MHdataframe[i,"Post pred SD"] ) )  , 3)


      MHdataframe[i,"FUT.A"]  = round( pnorm(1, mean = MHdataframe[i,"Post Mean"], sd =  MHdataframe[i,"Post pred SD"] ) - pnorm(-1, mean = MHdataframe[i,"Post Mean"], sd =  MHdataframe[i,"Post pred SD"] ) , 3)


      MHdataframe[i,"FUT.B+"] = round( ( pnorm( 1.5,  MHdataframe[i,"Post Mean"], MHdataframe[i,"Post pred SD"] ) - pnorm(1.0, MHdataframe[i,"Post Mean"], MHdataframe[i,"Post pred SD"] ) )  , 3)


      MHdataframe[i,"FUT.C+"] = round( ( pnorm( 100,  MHdataframe[i,"Post Mean"], MHdataframe[i,"Post pred SD"] ) - pnorm(1.5, MHdataframe[i,"Post Mean"], MHdataframe[i,"Post pred SD"] ) ) , 3)


    }

    return(MHdataframe)

  }

  #RUN THE EB FUNCTION

  MH_collection1 <- MHfunEB(MH_collection1)

  #######################################################################
  #PURIFICATION STEP: REMOVE ALL ITEMS THAT EXHIBIT C-LEVEL DIF AND RERUN
  #######################################################################

  if(purification == TRUE & any(( MH_collection1$"ETS Rule corr." == "C" )) ){

    print("Running purification")

    #IDENTIFY ITEMS WITH C-LEVEL DIF

    cleveldif <- which( MH_collection1$"ETS Rule corr." == "C"  )

    nocleveldif <- which( !(MH_collection1$"ETS Rule corr." == "C") )

    #REMOVE ITEMS THAT DISPLAY c-LEVEL DIF FROM ITEM FILE

    item.difficulty <- item.difficulty[-cleveldif]
    item.names <- item.names[-cleveldif]

    #RECALCULATE EXPECTED SCORE WITHOUT C-LEVEL ITEMS

    alldata[,"expected"] <- apply(sapply(item.difficulty, function(x) 1/(1+exp(-(person.ability-x)))),1,sum)

    #REMOVE ITEMS THAT DISPLAY C-LEVEL DIF FROM SCORED RESPONSE FILE

    alldata <- alldata[,-cleveldif]

    #NUMBER OF ITEMS DECREASES BY THE NUMBER OF C-LEVE DIF ITEMS DETECTED

    its <- its - length(cleveldif)

    #CREATE A NEW DATA FRAME FOR PURIFICATION

    MH_collection2 = data.frame(1:length(nocleveldif))

    #RUN THE MH FUNCTION

    MH_collection2 <- MHfun(MH_collection2)

    #RUN THE EB FUNCTION

    MH_collection2 <- MHfunEB(MH_collection2)

    #REPLACE C-LEVEL DIF ITEMS IN THE FINAL OUTPUT AFTER PURIFICATION OF NON-C-LEVEL DIF ITEMS

    MH_collection2 <- rbind(MH_collection2, MH_collection1[c(cleveldif),])

    # REORDER ITEMS TO ORIGINAL ORDER
    MH_collection2 <- MH_collection2[order(match(MH_collection2[,"Names"],MH_collection1[,"Names"])),]
    
    # FIX ITEM ORDER AGAIN
    MH_collection2[,"Item"] <- 1:(its+length(cleveldif))

    if( !is.null(outfilename) ){
      write.table(MH_collection2,outfilename,sep=",",row.names=FALSE);return(MH_collection2)
    }else{return(MH_collection2)}


  }else{
    if( !is.null(outfilename) ){
      write.table(MH_collection1,outfilename,sep=",",row.names=FALSE);return(MH_collection1)
    }else{return(MH_collection1) }


  }


}

