# An Application for studying Fentanyl clustering using Mass Spectra
# Arun Moorthy arun.moorthy@nist.gov // arunsmoorthy@gmail.com
# ==============================================================================

rm(list=ls())

# ==============================================================================
# Arun's Personal Functions and Routines
# ==============================================================================
source("Routines/asm-externalPackages.R")
source("Routines/asm-colours.R")

source("Functions/asm-PeakListCreator.R")
source("Functions/asm-TopPeaksSelector.R")
source("Functions/asm-SpectraPadder.R")
source("Functions/asm-SpectralSimilarity.R")
source("Functions/asm-simpleMF.R")
source("Functions/asm-hMF-IntensitySplitter.R")
source("Functions/asm-hybridMF.R")
source("Functions/asm-headtail_spectrumPlotter.R")
source("Functions/asm-simple-mw.R")
source("Functions/asm-strucPlotter-MODULAR-NEW.R")
source("Functions/asm-FentMassEstimator.R")
source("Functions/asm-LoadQuerySpectrum.R")
source("Functions/asm-MDS.R")
source("Functions/asm-MuIndex.R")
source("Functions/asm-cosineSim.R")


# ==============================================================================
# External Packages (Shiny)
#  =============================================================================
library(shiny)
library(shinythemes)
# ==============================================================================
library(MASS)

  # ==============================================================================
# Load Spectral Data (R Data Tables)
# ==============================================================================
stime <- system.time({
  FentanylLibrary<- readRDS("Data/SingleModFentanylLibrary.RDS")
})[3]
ReadInDTLibraryTime = stime;

# Load Structure Data (text file)
stime <- system.time({
  rawData = readLines("Data/asm-SingleModFentanylStructureDefinitions.txt")
  NameIndex = which(rawData=="$$$$")+1;
  Name = rawData[NameIndex]
  
  B_index = which(rawData=="[[B]]")
  A_index = which(rawData=="[[A]]")
  C_index = which(rawData=="[[C]]")
  D_index = which(rawData=="[[D]]")
  E_index = which(rawData=="[[E]]")
})[3]
ReadInDTLibraryTime = stime;

# Load Fent Map (R Data Structure)
FentMap = readRDS("Data/FentMap.RData")

# Default Values;
threshold_ab = 0;
ab_scaling = 0.5;
nSpectra = dim(FentanylLibrary)[1]
# ==============================================================================


# ==============================================================================
# Define UI for application 
# ==============================================================================
ui <- fixedPage(theme=shinytheme("flatly"),
                
      tags$head(tags$style(HTML('.modal.in .modal-dialog{
                                          width:100%;
                                          height:100%;
                                          margin:0px;
                                          }
                                          
                                          .modal-content{
                                          width:90%;
                                          height:90%;
                                          }
                                          '))),
      
  
      
  tabsetPanel(

    tabPanel("Search Results",
  
      br(),    
      titlePanel("Fentanyl Classifier"),
      # br(),
      div(p("The Fentanyl Classifier is a prototype implementation of \"augmented mass spectral library searching\". The software was designed for demonstration purposes. The authors cannot
      guarantee the accuracy of results generated using the Fentanyl Classifier, and cannot validate
      claims of others using this software."),
          style = "color:black;font-size:12pt"),
      
      # Input: Select a file 
      fileInput("file", "Choose Query Spectrum (MSP File)",
                multiple = FALSE,
                accept = c(".MSP")
      ),
      
      fixedRow(
          column(width=5,
            div(htmlOutput("QuerySentence"),style = "color:black;font-size:8pt"),
            plotOutput("QueryStruc", height = "280px", width="380px")
          ),
     
          column(width=7,
            htmlOutput("MWInfo",style="color:black;font-size:8pt"),    
            plotOutput("distPlot",height="300px",width="700px")
          )
        ),         

      
      fixedRow(   
          column(width=5,
            tabsetPanel(
                tabPanel("Hit List",
                   div(DT::dataTableOutput("HitList"), style = "font-size:70%")
                ),
                tabPanel("Hit Map",
                   plotOutput("MDSPlot",click="plot1_click",dblclick = "plot1_dblclick",
                              brush = brushOpts(id = "plot1_brush",resetOnNew = TRUE),
                              height="420px",width="420px"),
                   img(src='FentTypeLegend.png',align="center",width="420px")
                ), 
                tabPanel("Fentanyl Structure Definition", 
                         br(),
                         div(p("Fentanyl analogs can be described by the type and location of the structural 
                               modifications by which they differ from fentanyl. For example, alpha methyl 
                               fentanyl contains a methyl addition on the alpha position of modification site \"b\". 
                               The defined modification sites and structural scaffold were interpreted from the 
                               definitions provided in Fed Regist. 2018 Feb 6;83(25):5188-92."),style="text-align:justify;color:black;font-size:8pt"),
                         br(),
                         HTML('<p><img src="Pictures/Fentanyl.png", width="250px", align = "center"/></p>')
                )
            )
        
          ),
      
      
          column(width=2,
            br(),
            br(),
            br(),
            br(),br(),br(),
            htmlOutput("LibraryCompound",style="color:black;font-size:12pt"),
            div(tableOutput("LibraryData"), style = "color:black;font-size:8pt")
          ),
       
          column(width=5,
            plotOutput("LibraryStruc", height = "350px", width="550px")
          )        
      ),
      
   
   # Horizontal line 
   tags$hr()
   
   # tags$hr(),
   # fixedRow(
   #   column(width = 12,
   #          # Copyright messages
   #          div(p("(c) 2018",a("Mass Spectrometry Data Center", href="http://chemdata.nist.gov",target="_blank"), "at the National Institute of Standards and Technology."),
   #              br(), style = "color:black;font-size:6pt")
   #    )
   #   
   #  )
    
   ), 
  
  
  
  
  tabPanel("Instructions",
           br(),   
           titlePanel("Fentanyl Classifier (Instructions)"),
           # Horizontal line 
           tags$hr(),
           
           HTML('<p>1. Submit experimental spectrum as MSP file. </p>'),
           HTML('<p><img src="Pictures/Picture1.png", width="700px", align = "center"/></p>'),
           br(),
           HTML('<p>2. Review proposed structure generated through library search and automated interpretation. </p>'),
           HTML('<p><img src="Pictures/Picture2.png", width="700px", align = "center"/></p>'),
           br(),
           HTML('<p>3. Manually inspect library compounds, spectra, and search results.</p>'),
           HTML('<p><img src="Pictures/Picture3.png", width="700px", align = "center"/></p>'),
           
           # Horizontal line 
           tags$hr()
           # tags$hr(),
           # 
           # fixedRow(
           #   column(width = 12,
           #          # Copyright messages
           #          div(p("(c) 2018",a("Mass Spectrometry Data Center", href="http://chemdata.nist.gov",target="_blank"), "at the National Institute of Standards and Technology."),
           #              br(), style = "color:black;font-size:6pt"))
           # )
  )
)    
)



# ==============================================================================


# ==============================================================================
# Define server logic for application
# =============================================================================
server <- function(input, output,session) {
  
  values <- reactiveValues()
  
  values$mwisest <- 0; 
  values$mwused <- 0;
  values$querySpectrum <- 0;
  values$hitselection <- 1; 
  values$SRI <- 0;
  
  values$notFentanyl <- 1;
  values$isFentanyl <- 0;
  values$Type1 <- 0;
  values$Type2 <-0;


  trigger1 <- eventReactive(input$file$datapath,ignoreInit = F, {
     
     InpDirName = input$file$datapath  
     a = asm_load_QuerySpectrum(InpDirName,values);
     
     mwq = a[1];
     values$mwused <- mwq;
     
     q = a[2:length(a)];
     values$querySpectrum <- q;
     
     # Forward Library Search
     nSpectra = dim(FentanylLibrary)[1]
     Results = matrix(0L,nrow=nSpectra,ncol=3)
     
     withProgress(message = 'Forward Library Search', value = 0, {
     
     for(j in 1:nSpectra){
       
       l = FentanylLibrary[j,PeakList];
       mwl = FentanylLibrary[j,MW];
       
       Results[j,1:3] = asm_SpectralSimilarityMeasures(q,mwq,l,mwl,threshold_ab,ab_scaling,wplot=0)
      
       incProgress(j/nSpectra) 
     }

     })
     
     HitList = as.data.table(cbind(FentanylLibrary[,.(Name,ModType)],Results))
     names(HitList)[3] = "hMF";
     names(HitList)[4] = "sMF";
     names(HitList)[5] = "DeltaMass"
    
     return(HitList)
     
  })
  
  
  trigger2 <- eventReactive(input$file$datapath,ignoreInit = F, {

    mwq = values$mwused;
    q = values$querySpectrum; 
    
    
    # Backward Library Search
    nSpectra = dim(FentanylLibrary)[1]
    Results = matrix(0L,nrow=nSpectra,ncol=3)
    
    withProgress(message = 'Backward Library Search', value = 0, {
      
    for(j in 1:nSpectra){
      
      l = FentanylLibrary[j,PeakList];
      mwl = FentanylLibrary[j,MW];
      
      Results[j,1:3] = asm_SpectralSimilarityMeasures(l,mwl,q,mwq,threshold_ab,ab_scaling,wplot=0)
     
       incProgress(j/nSpectra) 
    }
      
    })
    
    HitList = as.data.table(cbind(FentanylLibrary[,.(Name,ModType)],Results))
    names(HitList)[3] = "hMF";
    names(HitList)[4] = "sMF";
    names(HitList)[5] = "DeltaMass"
    
    return(HitList)
    
  })
  
  
  ranges <- reactiveValues(tx = c(-1,1), ty = c(-1,1))
    
  
  output$MDSPlot <- renderPlot({
     
     
     HitListF = trigger1();
     HitListB = trigger2();
     
     result = asm_MDS(HitListF,HitListB,FentMap);
     x = result[,1];
     y = result[,2];
     
     IDs = match(FentanylLibrary[,Name],HitListF[order(-hMF),Name])
     
     cols = c(col2,col3,col4,col5,col6,col1,"white")     
     CID = c(7,match(FentanylLibrary[,ModType],LETTERS))
     sID = c(1e-8,(HitListF[,hMF]/999));
     aID = 0.25*sID;
     
     plot(x,y,
          pch=19,
          xlim=ranges$tx,ylim=ranges$ty,
          col=alpha(cols[CID],aID),
          cex=4*sID,
          xlab=" ",ylab=" ");
     text(x[-1],y[-1],IDs,cex=0.75*sID[-1])

     matpoints(x[1],y[1],pch=18,col=colGold,cex=3)
     
     
     fentID = c(1,14);
     p = x[-fentID];
     q = y[-fentID];
    
     forClust = cbind(p,q);
     kclust = kmeans(forClust,3,algorithm="Hartigan-Wong",nstart=20);
     
     index = which(kclust$cluster==1)
     dataEllipse(forClust[index,1],forClust[index,2],
                 levels=0.95, 
                 cex=0.1,
                 pch=3,
                 col="black",
                 lwd=0.5,
                 lty=3,
                 plot.points=FALSE,
                 add=TRUE)

     index = which(kclust$cluster==2)
     dataEllipse(forClust[index,1],forClust[index,2],
                 levels=0.95,
                 pch=3,
                 col="black",
                 lwd=0.5,
                 lty=3,
                 plot.points=FALSE,
                 add=TRUE)

     index = which(kclust$cluster==3)
     dataEllipse(forClust[index,1],forClust[index,2],
                 levels=0.95,
                 pch=3,
                 col="black",
                 lwd=0.5,
                 lty=3,
                 plot.points=FALSE,
                 add=TRUE)
     
    })
    
  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$plot1_dblclick, {
      brush <- input$plot1_brush
      if (!is.null(brush)) {
        ranges$tx <- c(brush$xmin, brush$xmax)
        ranges$ty <- c(brush$ymin, brush$ymax)
        
      } else {
        ranges$tx <- c(-1,1)
        ranges$ty <- c(-1,1)
      }
  })  
    
  
  observeEvent(input$HitList_rows_selected, {
    
    k = input$HitList_rows_selected 
    values$hitselection <- k
    
  })
  
  
  observeEvent(input$plot1_click, {
   
    HitListF = trigger1();
    HitListB = trigger2()
    
    result = asm_MDS(HitListF,HitListB,FentMap);
    x = result[,1];
    y = result[,2];
    
    ax = input$plot1_click$x
    ay = input$plot1_click$y
    
    dist = sqrt((x[-1]-ax)^2+(y[-1]-ay)^2);
  
    j = which.min(dist);
    
    data = HitListF[order(-hMF),Name]

    k = which(data == FentanylLibrary[j,Name])

    values$hitselection <- k
    
  })
  

  output$HitList <- DT::renderDataTable(DT::datatable({
     
     withProgress(message = 'Generating HitList', value = 0.5, {
     HitListF = trigger1();
     HitListB = trigger2()
     
     result = asm_MDS(HitListF,HitListB,FentMap);
     x = result[,1];
     y = result[,2];
     
     fentID = c(1,14);
     p = x[-fentID];
     q = y[-fentID];
     
     forClust = cbind(p,q);
     kclust = kmeans(forClust,3,algorithm="Hartigan-Wong",nstart=20);
     
     index1 = which(kclust$cluster==1);
     index2 = which(kclust$cluster==2);
     index3 = which(kclust$cluster==3);
     oclusts = numeric(42);
     oclusts[index1] = 1;
     oclusts[index2] = 2;
     oclusts[index3] = 3;
     
     clusts = c(oclusts[1:12],4,oclusts[13:43])
     
     withProgress(message = 'Computing Spectra Relatedness Indices', value = 0, {
     
     Confidence = numeric(nSpectra);
     for(i in 1:nSpectra){
       MFU = HitListF[i,hMF]/999;
       MDU = 1 - (sqrt((x[i+1]-x[1])^2 + (y[i+1]-y[1])^2)/sqrt(8))
       MDU = max(0,MDU)
       
       Confidence[i] =  round((MFU*MDU) ,2)
       
       incProgress(i/nSpectra)
     }

     })
     
     muIndex = asm_muIndex(HitListF,HitListB,FentMap);
     HitList = as.data.table(cbind(HitListF,Confidence,clusts))
     HitList = HitList[order(-hMF)]
     HitList = HitList[,.(Name,ModType,hMF,DeltaMass,Confidence,clusts,muIndex)];
     
     values$SRI <- Confidence;

     names(HitList)[2] = "Mod.";
     names(HitList)[4] = "DM" ;
     names(HitList)[5] = "SRI";
     names(HitList)[6] = "Cl.";
     names(HitList)[7] = "mu-Index"
     
    incProgress(0.5)
     
     })
     
     return(HitList[,c(1,3:5)])
     
  }),server = FALSE, selection=list(mode='single',selected=values$hitselection),
  options = list(lengthMenu = c(10, 20, 50), pageLength = 10))

   
  output$distPlot <- renderPlot({
     
     mwq = values$mwused; 
     q = values$querySpectrum; 
     
     HitList = trigger1();
     HitList = HitList[order(-hMF)]
     
     k = values$hitselection  
     
     i = which(FentanylLibrary[,Name]==HitList[k,Name])

     l = FentanylLibrary[i,PeakList];
     mwl = FentanylLibrary[i,MW];
     
     asm_SpectralSimilarityMeasures(q,mwq,l,mwl,threshold_ab,ab_scaling,wplot=1)
     
  })
   
   
  output$MWInfo <- renderText({
     trigger1();
     
     if(isolate(values$mwisest)==1){
       Message = paste("<p>No molecular weight information available from Query MSP.  An estimate of ", isolate(values$mwused), " Da was generated by adding 91 Da to the highest mass with intensity greater than 800.</a></p>",sep="");
     } else {
       Message = paste("<p>Molecular weight information was contained in the Query MSP. The nominal MW of ", isolate(values$mwused)," Da was employed in the Hybrid Search.</p>",sep="");
     }
     return(Message)
  })
  
   
  output$LibraryCompound <- renderText({
     trigger1()
     return("<b><u>Library Compound:</u></b>")
  })

  
  output$LibraryData <- renderTable({
     HitListF = trigger1();
     HitListB = trigger2()
     
     result = asm_MDS(HitListF,HitListB,FentMap);
     x = result[,1];
     y = result[,2];
     
     IDs = match(FentanylLibrary[,Name],HitListF[order(-hMF),Name])
     
     HitList = HitListF
     
     HitList = HitList[order(-hMF)]

     k = values$hitselection  
     
     i = which(FentanylLibrary[,Name]==HitList[k,Name])
     
     name = FentanylLibrary[i,Name];
     modType = FentanylLibrary[i,ModType];
     formula = FentanylLibrary[i,Formula]
     
     inchi = FentanylLibrary[i,InChIKey];
     mw = FentanylLibrary[i,MW];
     emass = FentanylLibrary[i,ExactMass];

     Message = as.data.frame(array(0,dim=c(5,1)));
     rownames(Message) <- c("Name:","Formula:","Exact Mass:","MW:","InChIKey:");
     Message[1,1] = name
     Message[2,1] = formula
     Message[3,1] = round(as.numeric(emass),3)
     Message[4,1] = mw
     Message[5,1] = paste("<a href='http://www.ncbi.nlm.nih.gov/sites/entrez?cmd=search&db=pccompound&term=\"",inchi,"\"[InChIKey]' target='_blank'>", inchi, "</a>")
     return(Message)
                
  },rownames=TRUE,colnames=FALSE,sanitize.text.function = function(x) x)
  
   
  output$LibraryStruc <- renderPlot({
     HitList = trigger1();
     HitList = HitList[order(-hMF)]

     i = values$hitselection
     k = which(FentanylLibrary[,Name]==HitList[i,Name])
     
     eval(parse(text=rawData[(B_index[k]+1):(B_index[k]+6)]))
     a_blockInfo = list(cbind(x,y,n));
     b_blockInfo = list(cbind(i,j,t));
     B = as.data.table(cbind(a_blockInfo,b_blockInfo));
     
     eval(parse(text=rawData[(A_index[k]+1):(A_index[k]+6)]))
     a_blockInfo = list(cbind(x,y,n));
     b_blockInfo = list(cbind(i,j,t));
     A = as.data.table(cbind(a_blockInfo,b_blockInfo));
     
     eval(parse(text=rawData[(C_index[k]+1):(C_index[k]+7)]))
     a_blockInfo = list(cbind(x,y,n));
     b_blockInfo = list(cbind(i,j,t,l));
     C = as.data.table(cbind(a_blockInfo,b_blockInfo));
     
     eval(parse(text=rawData[(D_index[k]+1):(D_index[k]+6)]))
     a_blockInfo = list(cbind(x,y,n));
     b_blockInfo = list(cbind(i,j,t));
     D = as.data.table(cbind(a_blockInfo,b_blockInfo));
     
     eval(parse(text=rawData[(E_index[k]+1):(E_index[k]+6)]))
     a_blockInfo = list(cbind(x,y,n));
     b_blockInfo = list(cbind(i,j,t));
     E = as.data.table(cbind(a_blockInfo,b_blockInfo));
     
     asm_strucPlotter_2(B,A,C,D,E)
     
  })
  
  
  output$QuerySentence <- renderText({
    trigger1();
    
    Message = "<p>Potential structure based on library search results. <b>Disclaimer: The authors do not guarantee <br>the accuracy of this result or claims of others based on results generated using this tool.</b></p>"
    
    return(Message)
  })
  

  output$QueryStruc <- renderPlot({
    HitListF = trigger1();
    HitListB = trigger2();
    
    result = asm_MDS(HitListF,HitListB,FentMap);
    x = result[,1];
    y = result[,2];
    
    QCenter = result[1,];
    
    fentID = c(1,14);
    p = x[-fentID];
    q = y[-fentID];
    
    forClust = cbind(p,q);
    kclust = kmeans(forClust,3,algorithm="Hartigan-Wong",nstart=20);
    
    index1 = which(kclust$cluster==1);
    index2 = which(kclust$cluster==2);
    index3 = which(kclust$cluster==3);
    oclusts = numeric(42);
    oclusts[index1] = 1;
    oclusts[index2] = 2;
    oclusts[index3] = 3;
    
    clusts = c(oclusts[1:12],4,oclusts[13:43])
    
    clustCenterG1 = kclust$centers[1,]
    clustCenterG2 = kclust$centers[2,]
    clustCenterG3 = kclust$centers[3,]
    
    clustSize1 = kclust$size[1];
    clustSize2 = kclust$size[2];
    clustSize3 = kclust$size[3];
    
    dG1center = sqrt((x[1]-clustCenterG1[1])^2 + (y[1]-clustCenterG1[2])^2);
    # print(dG1center)
    
    dG2center = sqrt((x[1]-clustCenterG2[1])^2 + (y[1]-clustCenterG2[2])^2);
    # print(dG2center)
    
    dG3center = sqrt((x[1]-clustCenterG3[1])^2 + (y[1]-clustCenterG3[2])^2);
    # print(dG3center)
    
    dCenters = c(dG1center,dG2center,dG3center);
    
    SRI <- values$SRI;

    HitListF = cbind(HitListF,SRI,clusts);
    HitList = HitListF[order(-SRI,-hMF)];
    relCompounds = HitList[(hMF>650)]; #only look at relevant compounds
    
    # print(relCompounds)
    
    # Are there any relevant compounds?
    nRelCompounds = dim(relCompounds)[1];
    if(nRelCompounds==0){
      plot(0,0,col="white",axes=FALSE,xlab="",ylab="")
      text(-0.5,0.25,"Not likely to be a fentanyl,")
      text(-0.52,0,   "Type I, or Type II analog.")
      return(0);
    }
    
    
    # Is the query likely an unmodified fentanyl?
    kF = which(relCompounds[,clusts]==4)
    nkF = length(kF); # 1 or 0
    if(nkF>0){
      if(relCompounds[kF,DeltaMass]==0 && relCompounds[kF,hMF]>850 && relCompounds[kF,SRI]>0.85){  # additional criteria to build confidence
        
        k = 13;
        
        eval(parse(text=rawData[(B_index[k]+1):(B_index[k]+6)]))
        a_blockInfo = list(cbind(x,y,n));
        b_blockInfo = list(cbind(i,j,t));
        B = as.data.table(cbind(a_blockInfo,b_blockInfo));
        
        eval(parse(text=rawData[(A_index[k]+1):(A_index[k]+6)]))
        a_blockInfo = list(cbind(x,y,n));
        b_blockInfo = list(cbind(i,j,t));
        A = as.data.table(cbind(a_blockInfo,b_blockInfo));
        
        eval(parse(text=rawData[(C_index[k]+1):(C_index[k]+7)]))
        a_blockInfo = list(cbind(x,y,n));
        b_blockInfo = list(cbind(i,j,t,l));
        C = as.data.table(cbind(a_blockInfo,b_blockInfo));
        
        eval(parse(text=rawData[(D_index[k]+1):(D_index[k]+6)]))
        a_blockInfo = list(cbind(x,y,n));
        b_blockInfo = list(cbind(i,j,t));
        D = as.data.table(cbind(a_blockInfo,b_blockInfo));
        
        eval(parse(text=rawData[(E_index[k]+1):(E_index[k]+6)]))
        a_blockInfo = list(cbind(x,y,n));
        b_blockInfo = list(cbind(i,j,t));
        E = as.data.table(cbind(a_blockInfo,b_blockInfo));
        
        asm_strucPlotter_2(B,A,C,D,E)
        
        return(0)
        
      }
    }
    

    # If nkF == 1 and the match factor is high, that means fentanyl is in the 
    # relevant compound list. In that case, either:
    # (1) it is a type 1 fentanyl analog so there will only be one modification site. Or
    # (2) it is a type 2 fentanyl analog with two modifications on site A and B
    
    
    if(nkF>0){
      
    ku = order(dCenters)[1]; #print(ku);
    
    k1 = which(relCompounds[,clusts]==ku);
    m1 = relCompounds[ku,ModType][1]; #print(m1)
    
      i = which(relCompounds[k1,DeltaMass]==0 && relCompounds[k1,hMF]>850 && relCompounds[k1,SRI]>0.85);
      
      if (length(i)>0){
        j = i[1];
        k = which(FentanylLibrary[,Name]==relCompounds[k1[j],Name])
        
        eval(parse(text=rawData[(B_index[k]+1):(B_index[k]+6)]))
        a_blockInfo = list(cbind(x,y,n));
        b_blockInfo = list(cbind(i,j,t));
        B = as.data.table(cbind(a_blockInfo,b_blockInfo));
        
        eval(parse(text=rawData[(A_index[k]+1):(A_index[k]+6)]))
        a_blockInfo = list(cbind(x,y,n));
        b_blockInfo = list(cbind(i,j,t));
        A = as.data.table(cbind(a_blockInfo,b_blockInfo));
        
        eval(parse(text=rawData[(C_index[k]+1):(C_index[k]+7)]))
        a_blockInfo = list(cbind(x,y,n));
        b_blockInfo = list(cbind(i,j,t,l));
        C = as.data.table(cbind(a_blockInfo,b_blockInfo));
        
        eval(parse(text=rawData[(D_index[k]+1):(D_index[k]+6)]))
        a_blockInfo = list(cbind(x,y,n));
        b_blockInfo = list(cbind(i,j,t));
        D = as.data.table(cbind(a_blockInfo,b_blockInfo));
        
        eval(parse(text=rawData[(E_index[k]+1):(E_index[k]+6)]))
        a_blockInfo = list(cbind(x,y,n));
        b_blockInfo = list(cbind(i,j,t));
        E = as.data.table(cbind(a_blockInfo,b_blockInfo));
        
        asm_strucPlotter_2(B,A,C,D,E)
        
        return(0)
        
      } else {
        
          k = 13
          
          eval(parse(text=rawData[(B_index[k]+1):(B_index[k]+6)]))
          a_blockInfo = list(cbind(x,y,n));
          b_blockInfo = list(cbind(i,j,t));
          B = as.data.table(cbind(a_blockInfo,b_blockInfo));
          
          eval(parse(text=rawData[(A_index[k]+1):(A_index[k]+6)]))
          a_blockInfo = list(cbind(x,y,n));
          b_blockInfo = list(cbind(i,j,t));
          A = as.data.table(cbind(a_blockInfo,b_blockInfo));
          
          eval(parse(text=rawData[(C_index[k]+1):(C_index[k]+7)]))
          a_blockInfo = list(cbind(x,y,n));
          b_blockInfo = list(cbind(i,j,t,l));
          C = as.data.table(cbind(a_blockInfo,b_blockInfo));
          
          eval(parse(text=rawData[(D_index[k]+1):(D_index[k]+6)]))
          a_blockInfo = list(cbind(x,y,n));
          b_blockInfo = list(cbind(i,j,t));
          D = as.data.table(cbind(a_blockInfo,b_blockInfo));
          
          eval(parse(text=rawData[(E_index[k]+1):(E_index[k]+6)]))
          a_blockInfo = list(cbind(x,y,n));
          b_blockInfo = list(cbind(i,j,t));
          E = as.data.table(cbind(a_blockInfo,b_blockInfo));
          
          if (m1 == "A"){
            k = 45;
            eval(parse(text=rawData[(A_index[k]+1):(A_index[k]+6)]))
            a_blockInfo = list(cbind(x,y,n));
            b_blockInfo = list(cbind(i,j,t));
            A = as.data.table(cbind(a_blockInfo,b_blockInfo));
          } else if (m1 == "B"){
            k=45;
            eval(parse(text=rawData[(B_index[k]+1):(B_index[k]+6)]))
            a_blockInfo = list(cbind(x,y,n));
            b_blockInfo = list(cbind(i,j,t));
            B = as.data.table(cbind(a_blockInfo,b_blockInfo));
          } else if (m1 == "C"){
            k=45;
            eval(parse(text=rawData[(C_index[k]+1):(C_index[k]+7)]))
            a_blockInfo = list(cbind(x,y,n));
            b_blockInfo = list(cbind(i,j,t,l));
            C = as.data.table(cbind(a_blockInfo,b_blockInfo));
          } else if (m1 == "D"){
            k=45;
            eval(parse(text=rawData[(D_index[k]+1):(D_index[k]+6)]))
            a_blockInfo = list(cbind(x,y,n));
            b_blockInfo = list(cbind(i,j,t));
            D = as.data.table(cbind(a_blockInfo,b_blockInfo));
          } else if (m1 == "E"){
            k=45;
            eval(parse(text=rawData[(E_index[k]+1):(E_index[k]+6)]))
            a_blockInfo = list(cbind(x,y,n));
            b_blockInfo = list(cbind(i,j,t));
            E = as.data.table(cbind(a_blockInfo,b_blockInfo));
          }
          
          asm_strucPlotter_2(B,A,C,D,E)
          
          return(0)
          
      }
    
    
        
    }
    
    
    # If nkF == 0, but nRelCompounds > 0, this means that this is likely a Type 2 (or greater) Fentanyl Analog and
    # so there will be two (or more) modification sites.
    
    a = unique(relCompounds[,clusts]);  # unique clusters in rel compounds list
    la = length(a)                      #number of unique clusters in rel compounds list
    
    # 
    if(la>1){
      k1 = which(relCompounds[,clusts]==order(dCenters)[1])[1];
      k2 = which(relCompounds[,clusts]==order(dCenters)[2])[1];
      
      m1 = relCompounds[k1,ModType];
      m2 = relCompounds[k2,ModType];
      
      k = 13
    
      eval(parse(text=rawData[(B_index[k]+1):(B_index[k]+6)]))
      a_blockInfo = list(cbind(x,y,n));
      b_blockInfo = list(cbind(i,j,t));
      B = as.data.table(cbind(a_blockInfo,b_blockInfo));
    
      eval(parse(text=rawData[(A_index[k]+1):(A_index[k]+6)]))
      a_blockInfo = list(cbind(x,y,n));
      b_blockInfo = list(cbind(i,j,t));
      A = as.data.table(cbind(a_blockInfo,b_blockInfo));
    
      eval(parse(text=rawData[(C_index[k]+1):(C_index[k]+7)]))
      a_blockInfo = list(cbind(x,y,n));
      b_blockInfo = list(cbind(i,j,t,l));
      C = as.data.table(cbind(a_blockInfo,b_blockInfo));
    
      eval(parse(text=rawData[(D_index[k]+1):(D_index[k]+6)]))
      a_blockInfo = list(cbind(x,y,n));
      b_blockInfo = list(cbind(i,j,t));
      D = as.data.table(cbind(a_blockInfo,b_blockInfo));
    
      eval(parse(text=rawData[(E_index[k]+1):(E_index[k]+6)]))
      a_blockInfo = list(cbind(x,y,n));
      b_blockInfo = list(cbind(i,j,t));
      E = as.data.table(cbind(a_blockInfo,b_blockInfo));
      
      
      k = which(FentanylLibrary[,Name]==relCompounds[k1,Name]);
      
      if (m1 == "A"){
        eval(parse(text=rawData[(A_index[k]+1):(A_index[k]+6)]))
        a_blockInfo = list(cbind(x,y,n));
        b_blockInfo = list(cbind(i,j,t));
        A = as.data.table(cbind(a_blockInfo,b_blockInfo));
      } else if (m1 == "B"){
        eval(parse(text=rawData[(B_index[k]+1):(B_index[k]+6)]))
        a_blockInfo = list(cbind(x,y,n));
        b_blockInfo = list(cbind(i,j,t));
        B = as.data.table(cbind(a_blockInfo,b_blockInfo));
      } else if (m1 == "C"){
        eval(parse(text=rawData[(C_index[k]+1):(C_index[k]+7)]))
        a_blockInfo = list(cbind(x,y,n));
        b_blockInfo = list(cbind(i,j,t,l));
        C = as.data.table(cbind(a_blockInfo,b_blockInfo));
      } else if (m1 == "D"){
        eval(parse(text=rawData[(D_index[k]+1):(D_index[k]+6)]))
        a_blockInfo = list(cbind(x,y,n));
        b_blockInfo = list(cbind(i,j,t));
        D = as.data.table(cbind(a_blockInfo,b_blockInfo));
      } else if (m1 == "E"){
        eval(parse(text=rawData[(E_index[k]+1):(E_index[k]+6)]))
        a_blockInfo = list(cbind(x,y,n));
        b_blockInfo = list(cbind(i,j,t));
        E = as.data.table(cbind(a_blockInfo,b_blockInfo));
      }
      
      k = which(FentanylLibrary[,Name]==relCompounds[k2,Name]);
      if (m2 == "A"){
        eval(parse(text=rawData[(A_index[k]+1):(A_index[k]+6)]))
        a_blockInfo = list(cbind(x,y,n));
        b_blockInfo = list(cbind(i,j,t));
        A = as.data.table(cbind(a_blockInfo,b_blockInfo));
      } else if (m2 == "B"){
        eval(parse(text=rawData[(B_index[k]+1):(B_index[k]+6)]))
        a_blockInfo = list(cbind(x,y,n));
        b_blockInfo = list(cbind(i,j,t));
        B = as.data.table(cbind(a_blockInfo,b_blockInfo));
      } else if (m2 == "C"){
        eval(parse(text=rawData[(C_index[k]+1):(C_index[k]+7)]))
        a_blockInfo = list(cbind(x,y,n));
        b_blockInfo = list(cbind(i,j,t,l));
        C = as.data.table(cbind(a_blockInfo,b_blockInfo));
      } else if (m2 == "D"){
        eval(parse(text=rawData[(D_index[k]+1):(D_index[k]+6)]))
        a_blockInfo = list(cbind(x,y,n));
        b_blockInfo = list(cbind(i,j,t));
        D = as.data.table(cbind(a_blockInfo,b_blockInfo));
      } else if (m2 == "E"){
        eval(parse(text=rawData[(E_index[k]+1):(E_index[k]+6)]))
        a_blockInfo = list(cbind(x,y,n));
        b_blockInfo = list(cbind(i,j,t));
        E = as.data.table(cbind(a_blockInfo,b_blockInfo));
      }
      
      
      asm_strucPlotter_2(B,A,C,D,E);
      
      return(0);

    } else if (la==1) {
      
      k1 = which(relCompounds[,clusts]==order(dCenters)[1])[1];
      m1 = relCompounds[k1,ModType];
      
      k2 = which(HitList[,clusts]==order(dCenters)[2])[1];
      m2 = HitList[k2,ModType];

      k = 13
      
      eval(parse(text=rawData[(B_index[k]+1):(B_index[k]+6)]))
      a_blockInfo = list(cbind(x,y,n));
      b_blockInfo = list(cbind(i,j,t));
      B = as.data.table(cbind(a_blockInfo,b_blockInfo));
      
      eval(parse(text=rawData[(A_index[k]+1):(A_index[k]+6)]))
      a_blockInfo = list(cbind(x,y,n));
      b_blockInfo = list(cbind(i,j,t));
      A = as.data.table(cbind(a_blockInfo,b_blockInfo));
      
      eval(parse(text=rawData[(C_index[k]+1):(C_index[k]+7)]))
      a_blockInfo = list(cbind(x,y,n));
      b_blockInfo = list(cbind(i,j,t,l));
      C = as.data.table(cbind(a_blockInfo,b_blockInfo));
      
      eval(parse(text=rawData[(D_index[k]+1):(D_index[k]+6)]))
      a_blockInfo = list(cbind(x,y,n));
      b_blockInfo = list(cbind(i,j,t));
      D = as.data.table(cbind(a_blockInfo,b_blockInfo));
      
      eval(parse(text=rawData[(E_index[k]+1):(E_index[k]+6)]))
      a_blockInfo = list(cbind(x,y,n));
      b_blockInfo = list(cbind(i,j,t));
      E = as.data.table(cbind(a_blockInfo,b_blockInfo));
      
      
      k = which(FentanylLibrary[,Name]==relCompounds[k1,Name]);
      
      if (m1 == "A"){
        eval(parse(text=rawData[(A_index[k]+1):(A_index[k]+6)]))
        a_blockInfo = list(cbind(x,y,n));
        b_blockInfo = list(cbind(i,j,t));
        A = as.data.table(cbind(a_blockInfo,b_blockInfo));
      } else if (m1 == "B"){
        eval(parse(text=rawData[(B_index[k]+1):(B_index[k]+6)]))
        a_blockInfo = list(cbind(x,y,n));
        b_blockInfo = list(cbind(i,j,t));
        B = as.data.table(cbind(a_blockInfo,b_blockInfo));
      } else if (m1 == "C"){
        eval(parse(text=rawData[(C_index[k]+1):(C_index[k]+7)]))
        a_blockInfo = list(cbind(x,y,n));
        b_blockInfo = list(cbind(i,j,t,l));
        C = as.data.table(cbind(a_blockInfo,b_blockInfo));
      } else if (m1 == "D"){
        eval(parse(text=rawData[(D_index[k]+1):(D_index[k]+6)]))
        a_blockInfo = list(cbind(x,y,n));
        b_blockInfo = list(cbind(i,j,t));
        D = as.data.table(cbind(a_blockInfo,b_blockInfo));
      } else if (m1 == "E"){
        eval(parse(text=rawData[(E_index[k]+1):(E_index[k]+6)]))
        a_blockInfo = list(cbind(x,y,n));
        b_blockInfo = list(cbind(i,j,t));
        E = as.data.table(cbind(a_blockInfo,b_blockInfo));
      }
      
      k = 45;
      if (m2 == "A"){
        eval(parse(text=rawData[(A_index[k]+1):(A_index[k]+6)]))
        a_blockInfo = list(cbind(x,y,n));
        b_blockInfo = list(cbind(i,j,t));
        A = as.data.table(cbind(a_blockInfo,b_blockInfo));
      } else if (m2 == "B"){
        eval(parse(text=rawData[(B_index[k]+1):(B_index[k]+6)]))
        a_blockInfo = list(cbind(x,y,n));
        b_blockInfo = list(cbind(i,j,t));
        B = as.data.table(cbind(a_blockInfo,b_blockInfo));
      } else if (m2 == "C"){
        eval(parse(text=rawData[(C_index[k]+1):(C_index[k]+7)]))
        a_blockInfo = list(cbind(x,y,n));
        b_blockInfo = list(cbind(i,j,t,l));
        C = as.data.table(cbind(a_blockInfo,b_blockInfo));
      } else if (m2 == "D"){
        eval(parse(text=rawData[(D_index[k]+1):(D_index[k]+6)]))
        a_blockInfo = list(cbind(x,y,n));
        b_blockInfo = list(cbind(i,j,t));
        D = as.data.table(cbind(a_blockInfo,b_blockInfo));
      } else if (m2 == "E"){
        eval(parse(text=rawData[(E_index[k]+1):(E_index[k]+6)]))
        a_blockInfo = list(cbind(x,y,n));
        b_blockInfo = list(cbind(i,j,t));
        E = as.data.table(cbind(a_blockInfo,b_blockInfo));
      }
      
      
      asm_strucPlotter_2(B,A,C,D,E);
      
      return(0);
      
    }

  })
  
  
  output$QueryData <- renderTable({
    
    HitListF = trigger1();
    HitListB = trigger2();
    
    result = asm_MDS(HitListF,HitListB,FentMap);
    x = result[,1];
    y = result[,2];
    
    QCenter = result[1,];
    
    fentID = c(1,14);
    p = x[-fentID];
    q = y[-fentID];
    
    forClust = cbind(p,q);
    kclust = kmeans(forClust,3,algorithm="Hartigan-Wong",nstart=20);
    
    index1 = which(kclust$cluster==1);
    index2 = which(kclust$cluster==2);
    index3 = which(kclust$cluster==3);
    oclusts = numeric(42);
    oclusts[index1] = 1;
    oclusts[index2] = 2;
    oclusts[index3] = 3;
    
    clusts = c(oclusts[1:12],4,oclusts[13:42])
    
    clustCenterG1 = kclust$centers[1,]
    clustCenterG2 = kclust$centers[2,]
    clustCenterG3 = kclust$centers[3,]
    
    clustSize1 = kclust$size[1];
    clustSize2 = kclust$size[2];
    clustSize3 = kclust$size[3];
    
    
    SRI <- values$SRI;
    
    HitListF = cbind(HitListF,SRI,clusts);
    HitList = HitListF[order(-hMF)];
    
    
    Message = as.data.frame(array(0,dim=c(5,1)));
    Message[1,1] = HitList[1,Name];
    Message[2,1] = HitList[1,ModType];
    Message[3,1] = HitList[1,SRI];
    Message[4,1] = HitList[1,Name];
    Message[5,1] = HitList[1,Name];
    
    return(Message)

  },rownames=TRUE,colnames=FALSE,sanitize.text.function = function(x) x)
  
  
  session$onSessionEnded(function() {
    stopApp()
  })


} # END SERVER LOGIC


# Run the application 
# ==============================================================================


# ==============================================================================
# Run Application
# ==============================================================================
shinyApp(ui = ui, server = server)
# ==============================================================================