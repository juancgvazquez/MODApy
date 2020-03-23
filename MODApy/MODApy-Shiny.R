#Req Libraries
library(shiny)
library(ConfigParser)
library(reticulate)
library(DT)
library(openxlsx)
library(shinycssloaders)
#python config -------------------------------------------------------------------
modapydir ='/home/darrow/miniconda3/envs/MODApy/lib/python3.7/site-packages/MODApy/'
cfgpath = paste0(modapydir,'config.ini')
logfile = paste0(modapydir,'logs/currentrun.log')
pipelog = paste0(modapydir,'logs/pipe_run.log')
pipeflag = paste0(modapydir,'logs/pipe.flag')
dlog = paste0(modapydir,'logs/downloads.log')
cfg = read.ini(cfgpath)
use_virtualenv("/home/darrow/miniconda3/envs/MODApy")
use_python("/home/darrow/miniconda3/envs/MODApy/bin/python3")
py_config()
MODApy<-import('MODApy')

# combo box options -------------------------------------------------------------------
patientsvcf <- result<-gsub('\\..*','',basename(list.files(cfg$PATHS$patientpath,pattern="\\.final.vcf",recursive = TRUE)))
patientsbam <- result<-gsub('\\..*','',basename(list.files(cfg$PATHS$patientpath,pattern="\\.bam",recursive = TRUE)))
patientsfastq <- result<-gsub('\\..*','',basename(list.files(cfg$PATHS$patientpath,pattern="\\.fastq.gz|\\.fastq|\\.fq|\\.fq.gz",recursive = TRUE)))
chromdbs <- gsub('\\..*','',basename(list.files(dirname(cfg$PATHS$dbpath),pattern="\\.csv",recursive = TRUE)))
updatepanels <<- function(){
  result<-gsub('.xlsx','',list.files(cfg$PATHS$panelspath,pattern='\\.xlsx$'))
  return(result)
}
panels <- gsub('.xlsx','',list.files(cfg$PATHS$panelspath,pattern='\\.xlsx$'))
# utils ---------------------------------------------------------------------------
getcommand <- function(input){
  cmd = ''
  switch (input$tabset,
          Pipeline={
            if(file.exists(paste(cfg$PATHS$patientpath, input$PatientPipe, "/",input$PatientPipe,"_1", ".fastq", sep=""))){
              cmd =  paste("pipeline -Pipeline",
                           input$Pipeline, "-FQ",
                           paste(input$PatientPipe, "/",input$PatientPipe,"_1", ".fastq", sep=""), "-FQ",
                           paste(input$PatientPipe, "/",input$PatientPipe,"_2", ".fastq", sep=""))
            }
            else if(file.exists(paste(cfg$PATHS$patientpath, input$PatientPipe, "/",input$PatientPipe,"_1", ".fastq.gz", sep=""))){
              cmd = paste("pipeline -Pipeline",
                          input$Pipeline, "-FQ",
                          paste(input$PatientPipe, "/",input$PatientPipe,"_1", ".fastq.gz", sep=""), "-FQ",
                          paste(input$PatientPipe, "/",input$PatientPipe,"_2", ".fastq.gz", sep=""))
            }
            else {
              cat("No fastq or gzipped fastq files found for that Patient",file=logfile,sep='\n')
            }
          },
          Single={
            cmd = paste("single -Panel", input$Panel, "-Patient", paste(input$PatientPanel, "/",input$PatientPanel,".final.vcf", sep=""))
          },
          Duos={
            cmd = paste("duos -Patient1", paste(input$Patient1D, "/",input$Patient1D,".final.vcf", sep=""), "-Patient2", paste(input$Patient2D, "/",input$Patient2D,".final.vcf", sep=""), '--Filter')
            venn=''
            if(input$vennplaceD == input$Patient1D){venn = "A"}
            else if(input$vennplaceD == input$Patient2D){venn = "B"}
            else if(input$vennplaceD == paste(input$Patient1D,input$Patient2D,sep=":")){venn = "A:B"}
            if(venn!=''){
              cmd = paste(cmd,"--VennPlace", venn)
            }
            if(input$PanelD != 'NONE'){
              cmd = paste(cmd, '--Panel', input$PanelD)
            }
          },
          Trios={
            venn=''
            cmd = cmd = paste("trios -Patient1", paste(input$Patient1T, "/",input$Patient1T,".final.vcf", sep=""), "-Patient2", paste(input$Patient2T, "/",input$Patient2T,".final.vcf", sep=""), "-Patient3", paste(input$Patient3T, "/",input$Patient3T,".final.vcf", sep=""), '--Filter')
            if(input$vennplaceT == input$Patient1T){venn = "A"}
            else if(input$vennplaceT == input$Patient2T){venn = "B"}
            else if(input$vennplaceT == input$Patient3T){venn = "C"}
            else if(input$vennplaceT == paste(input$Patient1T,input$Patient2T,sep=":")){venn = "A:B"}
            else if(input$vennplaceT == paste(input$Patient1T,input$Patient3T,sep=":")){venn = "A:C"}
            else if(input$vennplaceT == paste(input$Patient2T,input$Patient3T,sep=":")){venn = "B:C"}
            else if(input$vennplaceT == paste(input$Patient1T,input$Patient2T,input$Patient3T,sep=":")){venn = "A:B:C"}
            if(venn!=''){
              cmd = paste(cmd,"--VennPlace", venn)
            }
            if(input$PanelT != 'NONE'){
              cmd = paste(cmd, '--Panel', input$PanelT)
            }
          }
  )
  cmd <- gsub("[(]","\\\\(",cmd)
  cmd <- gsub("[)]","\\\\)",cmd)
  system2("MODApy", cmd ,wait=FALSE,stdout = FALSE,stderr = FALSE)
}

# ui -------------------------------------------------------------------
ui <- tagList(shinyjs::useShinyjs(),
              navbarPage("MODApy",
                         tabPanel("Analisis",
                                  sidebarLayout(
                                    sidebarPanel(
                                      width = 4,
                                      tabsetPanel(id="tabset", selected = "Single",
                                                  tabPanel("Single",
                                                           fluidRow(style='padding-left:20px;',
                                                                    br(),
                                                                    div(style="display: inline-block;vertical-align:center; width: 150px;",h4('Patient:')),
                                                                    div(style="display: inline-block;vertical-align:top; width: 300px;",selectInput(inputId = "PatientPanel", label = NULL, choices = patientsvcf)),
                                                                    div(style="display: inline-block;vertical-align:top; width: 50px;",actionButton('newpatient',label=NULL,icon=icon('plus')))
                                                           ),
                                                           fluidRow(style='padding-left:20px;',
                                                                    br(),
                                                                    div(style="display: inline-block;vertical-align:center; width: 150px;",h4('Panel:')),
                                                                    div(style="display: inline-block;vertical-align:top; width: 300px;",selectInput(inputId = "Panel", label = NULL, choices = panels)),
                                                                    div(style="display: inline-block;vertical-align:top; width: 50px;",actionButton('newpanel',label=NULL,icon=icon('plus'))))
                                                  ),
                                                  tabPanel("Duos",
                                                           fluidRow(style='padding-left:20px;',
                                                                    br(),
                                                                    div(style="display: inline-block;vertical-align:center; width: 150px;",h4('Patient 1:')),
                                                                    div(style="display: inline-block;vertical-align:top; width: 300px;",selectInput(inputId = "Patient1D", label = NULL, choices = patientsvcf)),
                                                                    div(style="display: inline-block;vertical-align:top; width: 50px;",actionButton('newpatient1D',label=NULL,icon=icon('plus')))
                                                           ),
                                                           fluidRow(style='padding-left:20px;',
                                                                    br(),
                                                                    div(style="display: inline-block;vertical-align:center; width: 150px;",h4('Patient 2:')),
                                                                    div(style="display: inline-block;vertical-align:top; width: 300px;",selectInput(inputId = "Patient2D", label = NULL, choices = patientsvcf)),
                                                                    div(style="display: inline-block;vertical-align:top; width: 50px;",actionButton('newpatient2D',label=NULL,icon=icon('plus')))
                                                           ),
                                                           fluidRow(style='padding-left:20px;',
                                                                    br(),
                                                                    div(style="display: inline-block;vertical-align:center; width: 150px;",h4('Panel:')),
                                                                    div(style="display: inline-block;vertical-align:top; width: 300px;",selectInput(inputId = "PanelD", label = NULL, choices = c('NONE',panels))),
                                                                    div(style="display: inline-block;vertical-align:top; width: 50px;",actionButton('newpanelD',label=NULL,icon=icon('plus')))
                                                           ),
                                                           fluidRow(style='padding-left:20px;',
                                                                    br(),
                                                                    div(style="display: inline-block;vertical-align:top; width: 150px;",h4('Venn Place:')),
                                                                    div(style="display: inline-block;vertical-align:top; width: 300px;",uiOutput("vennDuos"))
                                                           )
                                                  ),
                                                  tabPanel("Trios",
                                                           fluidRow(style='padding-left:20px;',
                                                                    br(),
                                                                    div(style="display: inline-block;vertical-align:center; width: 150px;",h4('Patient 1:')),
                                                                    div(style="display: inline-block;vertical-align:top; width: 300px;",selectInput(inputId = "Patient1T", label = NULL, choices = patientsvcf)),
                                                                    div(style="display: inline-block;vertical-align:top; width: 50px;",actionButton('newpatient1T',label=NULL,icon=icon('plus')))
                                                           ),
                                                           fluidRow(style='padding-left:20px;',
                                                                    br(),
                                                                    div(style="display: inline-block;vertical-align:center; width: 150px;",h4('Patient 2:')),
                                                                    div(style="display: inline-block;vertical-align:top; width: 300px;",selectInput(inputId = "Patient2T", label = NULL, choices = patientsvcf)),
                                                                    div(style="display: inline-block;vertical-align:top; width: 50px;",actionButton('newpatient2T',label=NULL,icon=icon('plus')))
                                                           ),
                                                           fluidRow(style='padding-left:20px;',
                                                                    br(),
                                                                    div(style="display: inline-block;vertical-align:center; width: 150px;",h4('Patient 3:')),
                                                                    div(style="display: inline-block;vertical-align:top; width: 300px;",selectInput(inputId = "Patient3T", label = NULL, choices = patientsvcf)),
                                                                    div(style="display: inline-block;vertical-align:top; width: 50px;",actionButton('newpatient3T',label=NULL,icon=icon('plus')))
                                                           ),
                                                           fluidRow(style='padding-left:20px;',
                                                                    br(),
                                                                    div(style="display: inline-block;vertical-align:center; width: 150px;",h4('Panel:')),
                                                                    div(style="display: inline-block;vertical-align:top; width: 300px;",selectInput(inputId = "PanelT", label = NULL, choices = c('NONE',panels))),
                                                                    div(style="display: inline-block;vertical-align:top; width: 50px;",actionButton('newpanelT',label=NULL,icon=icon('plus')))
                                                           ),
                                                           fluidRow(style='padding-left:20px;',
                                                                    br(),
                                                                    div(style="display: inline-block;vertical-align:top; width: 150px;",h4('Venn Place:')),
                                                                    div(style="display: inline-block;vertical-align:top; width: 300px;",uiOutput("vennTrios"))
                                                           )
                                                  )
                                      ),
                                      actionButton("buttonrun","Run"),
                                      actionButton("buttonlastcmd","Last Run Status"),
                                      shinyjs::disabled(downloadButton('downloadData','Download Result'))
                                    ),
                                    mainPanel(#style="background-color:blue",
                                      width = 6,
                                      h1('Command Output',style = "font-family: 'Courgette', cursive;
                                         font-weight: 500; line-height: 1.1; 
                                         color: #4d3a7d;"),
                                      htmlOutput("consoleout")
                                      )
                                    )),
                         tabPanel('VariantsDB',
                                  actionButton("buildDBbtn","Build DataBase"),
                                  div(style="display: inline-block;vertical-align:top; width: 300px;",selectInput(inputId = "ChromSel", label = NULL, choices = chromdbs)),
                                  actionButton("openDBbtn","Open Database"),
                                  actionButton("annotateFile",'Annotate File'),
                                  htmlOutput("dbout") %>% withSpinner(color="#0dc5c1"),
                                  DT::dataTableOutput("mytable")
                         ),
                         tabPanel("WES Pipelines",
                                  sidebarPanel(width=3,
                                    fluidRow(
                                      radioButtons('origin','Select Input File Type',choices = c('Fast Q'='fromfq','Processed BAM File'='frombam','Raw VCF'='fromvcf')
                                    ),
                                    conditionalPanel(
                                      'input.origin == "fromfq"',
                                      div(style="display: inline-block;vertical-align:top; width: 300px;",selectInput(inputId = "fqfile", label = NULL, choices = patientsfastq)),
                                      fileInput('fquploaded','Choose File to Upload',accept=c('fq.gzip','fastq.gzip','.fastq','.fq')),
                                      actionButton("runPipelinefq",'Run Pipeline')
                                    ),
                                    conditionalPanel(
                                      'input.origin == "frombam"',
                                      div(style="display: inline-block;vertical-align:top; width: 300px;",selectInput(inputId = "bamfile", label = NULL, choices = patientsbam)),
                                      fileInput('bamuploaded','Choose File to Upload',accept=c('.bam')),
                                      actionButton("runPipelinebam",'Run Pipeline')
                                    ),
                                    conditionalPanel(
                                      'input.origin == "fromvcf"',
                                      div(style="display: inline-block;vertical-align:top; width: 300px;",selectInput(inputId = "vcffile", label = NULL, choices = patientsvcf)),
                                      fileInput('vcfuploaded','Choose File to Upload',accept=c('.vcf')),
                                      actionButton("runPipelinevcf",'Run Pipeline')
                                      )
                                  )),
                                  mainPanel(
                                    conditionalPanel(
                                      'input.origin=="fromfq"',
                                      sidebarPanel(
                                        h2("Steps"),
                                        checkboxInput('trim','Trim with TrimGalore',value=TRUE),
                                        #checkboxInput('align','Align with BWA',value=TRUE),
                                        #checkboxInput('sort','Sort Generated Sam File',value=TRUE),
                                        checkboxInput('dedup','Mark Duplicates with Picard',value=TRUE),
                                        checkboxInput('recal','Recalibrate Bam File',value=TRUE),
                                        checkboxInput('callvars','Call Variants',value=TRUE),
                                        checkboxInput('annDBSNP','Annotate with DBSNP',value=TRUE),
                                        checkboxInput('ann1000GP3','Annotate with 1000GP3',value=TRUE),
                                        checkboxInput('annCLINVAR','Annotate with CLINVAR',value=TRUE),
                                        checkboxInput('annESP6500','Annotate with ESP6500',value=TRUE),
                                        checkboxInput('annGNOMAD','Annotate with GenomeAD',value=TRUE),
                                        checkboxInput('annDBNSFP','Annotate with DBNSFP',value=TRUE)
                                      ),
                                    mainPanel(
                                        h2("Full Pipeline Diagram"),
                                        img(src='wespipelinefull.png',width=700,align='center')
                                      )),
                                  conditionalPanel(
                                    'input.origin=="frombam"',
                                    sidebarPanel(
                                      h2("Steps"),
                                      checkboxInput('callvars','Call Variants',value=TRUE),
                                      checkboxInput('annDBSNP','Annotate with DBSNP',value=TRUE),
                                      checkboxInput('ann1000GP3','Annotate with 1000GP3',value=TRUE),
                                      checkboxInput('annCLINVAR','Annotate with CLINVAR',value=TRUE),
                                      checkboxInput('annESP6500','Annotate with ESP6500',value=TRUE),
                                      checkboxInput('annGNOMAD','Annotate with GenomeAD',value=TRUE),
                                      checkboxInput('annDBNSFP','Annotate with DBNSFP',value=TRUE)
                                    ),
                                    mainPanel(
                                      h2("Full Pipeline Diagram"),
                                      img(src='wespipelinebam.png',width=700,align='center')
                                    )
                                    ),
                                  conditionalPanel(
                                    'input.origin=="fromvcf"',
                                    sidebarPanel(
                                      h2("Steps"),
                                      checkboxInput('annDBSNP','Annotate with DBSNP',value=TRUE),
                                      checkboxInput('ann1000GP3','Annotate with 1000GP3',value=TRUE),
                                      checkboxInput('annCLINVAR','Annotate with CLINVAR',value=TRUE),
                                      checkboxInput('annESP6500','Annotate with ESP6500',value=TRUE),
                                      checkboxInput('annGNOMAD','Annotate with GenomeAD',value=TRUE),
                                      checkboxInput('annDBNSFP','Annotate with DBNSFP',value=TRUE)
                                    ),
                                    mainPanel(
                                      h2("Full Pipeline Diagram"),
                                      img(src='wespipelinevcf.png',width=700,align='center')
                                    )))
                                  #selectInput(inputId = "Pipeline", label = "Pipelines", choices = list.files(path=cfg$PATHS$pipelinespath)),
                                  #selectInput(inputId = "PatientPipe", label = "Patient", choices = list.dirs(
                                   # path = cfg$PATHS$patientpath ,full.names = FALSE,recursive = FALSE))
                         ),
                         tabPanel('About',
                                  h2('MODApy'),
                                  h3('Multi-Omics Data Analysis in Python - Shiny Frontend'),
                                  br(),
                                  p('MODApy is currently in development. Backend is developed in Python'),
                                  p('Frontend is developed in R through Shiny.'),
                                  p('Authors: Juan Carlos Vázquez - Elmer A. Fernández'),
                                  p('Bioscience Data Mining Group - Universidad Católica de Córdoba'),
                                  p('Centro de Investigación y Desarrollo en Inmunología y Enfermedades Infecciosas - CONICET'),
                                  actionButton("cfgButton",label=NULL,icon=icon('cog'))
                         )


              ),
              tags$head(tags$style(".modal-dialog{width:1000px}"))

)
# server -------------------------------------------------------------------
server <- function(input,output, session){
  dbfile = paste(dirname(cfg$PATHS$dbpath),'/',chromdbs[1],'.csv',sep="")
  if(file.exists(dbfile)){
    df1 <- read.csv(dbfile)
    output$mytable = DT::renderDataTable({df1},options=list(scrollX=TRUE))
  }
  rv <- reactiveValues(textstream = c(""), timer = reactiveTimer(1000),started=FALSE)
  rv2 <- reactiveValues(textstream2 = c(""), timer = reactiveTimer(1000),started=FALSE)
  rv3 <- reactiveValues(textstream3 = c(""), timer = reactiveTimer(1000),started=FALSE)
  downpath <- reactiveValues()
  newpatientModal <- function(failed = FALSE) {
    modalDialog(
      title = "Add New Patient",
      p('This will add a new patient, downlading the data from a url or a xlsx/xls file.
        Please either write the url or upload a file and press Submit. If both values are filled, will only download urls from files. 
        It will generate the folder under the Patient folder and download the .tar file'),
      textInput('url',NULL,value="",placeholder='Enter URL'),
      fileInput('file1','Choose File to Upload',accept=c('.xls','.xlsx')),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("addbtn","Add Patient")
      ))
  }
  annotateModal <- function(failed = FALSE) {
    modalDialog(
      title = "Annotate a file with Variants DB Freqs",
      p('This will annotate an excel file with the variant Frequencies available in Variants DB, both allelic and total'),
      fileInput('file1','Choose File to Annotate',accept=c('.xls','.xlsx')),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("annotatebtn","Annotate")
      ))
  }
  newpanelModal <- function(failed = FALSE) {
    modalDialog(
      title = "Add New Panel",
      p('This will add a new Panel. A panel consists in a list of genes'),
      p('Panel must be an Xlsx file with the same format as the example'),
      p(url <- a("Example Panel", href="data/Panel_Example.xlsx")),
      p('It is required to have a Sheet named \"GeneList\" with a column named \"GeneSymbol\ in it."'),
      p('Example can be downloaded and customized or upload a new file as long as the format is exactly equal.'),
      fileInput('panelupload','Choose Panel File to Upload',accept=c(".xlsx","xlsx","application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("addpanelbtn","Add Panel")
      ))
  }
  configModal <- function(failed = FALSE) {
    modalDialog(
      title = "Config",
      textInput("cfgcores","Number of Cores:",cfg$GENERAL$cores),
      textInput("cfgpatientspath","Patients Files Path:",cfg$PATHS$patientpath),
      textInput("cfgpanelspath","Paneles Files Path:",cfg$PATHS$panelspath),
      textInput("cfgresultspath","Results Files Path:",cfg$PATHS$resultspath),
      actionButton('savecfg',"Save Changes")
    )
  }
  pipelineRunningModal <- function(failed = FALSE) {
    modalDialog(
      title = "Warning",
      p('There is a pipeline already running'),
      p('Currently at:'),
      HTML(paste(readLines(pipelog),collapse='<br>'))
    )
  }
  observeEvent(c(input$newpatient,
                 input$newpatient1D,
                 input$newpatient2D,
                 input$newpatient1T,
                 input$newpatient2T,
                 input$newpatient3T), {
                   validate(need(input$newpatient > 0 | input$newpatient1D > 0 | input$newpatient2D > 0 |
                                   input$newpatient1T > 0 | input$newpatient2T > 0 | input$newpatient3T > 0, ''))
                   showModal(newpatientModal())
                 })
  observeEvent(c(input$newpanel,
                 input$newpanelD,
                 input$newpanelT), {
                   validate(need(input$newpanel > 0 | input$newpanelD > 0 | input$newpanelT > 0, ''))
                   showModal(newpanelModal())
                 })
  observeEvent(input$addpanelbtn,{
    if(tolower(tools::file_ext(input$panelupload$datapath)) == "xlsx"){
      file.copy(input$panelupload$datapath, paste(cfg$PATHS$panelspath,input$panelupload$name))
      removeModal()
      updateSelectInput(session,'Panel',choices=updatepanels())
      updateSelectInput(session,'PanelD',selected = 'NONE',choices=c('NONE',updatepanels()))
      updateSelectInput(session,'PanelT',selected = 'NONE',choices=c('NONE',updatepanels()))
    }
    else{
      removeModal()
      showModal(modalDialog(title='Error', 'Uploaded File must be an Excel file with xlsx extension'))
    }
  })
  observeEvent(input$buildDBbtn, {
    rv$textstream = ""
    rv$started<-TRUE
    system2('MODApy',args = 'variantsDB -buildDB',wait = FALSE,stdout = FALSE,stderr = FALSE)
  })
  observeEvent(input$openDBbtn, {
    dbfile = paste(dirname(cfg$PATHS$dbpath),'/',input$ChromSel,'.csv',sep="")
    if(file.exists(dbfile)){
      rv$textstream = ""
      rv$started<-FALSE
      df1 <- read.csv(dbfile)
      output$mytable = DT::renderDataTable({df1},options=list(scrollX=TRUE))
    }
    else{
      rv$textstream = "Variants file not found. Try to build variantsDB first."
    }

  })
  observeEvent(input$annotateFile, {
    showModal(annotateModal())
    })

  observeEvent(input$annotatebtn, {
    if(is.null(input$file1)){
      removeModal()
      modalDialog('No input data to download.')
    }
    else if(!(is.null(input$file1))){
      file.copy(input$file1$datapath, paste0("./", input$file1$name))
      withProgress(message='Annotating',value=0, {
        system2('MODApy',args = paste('variantsDB -annotate',paste0('./',input$file1$name)),wait = TRUE,stdout = FALSE,stderr = FALSE)%>% withSpinner(color="#0dc5c1")
        incProgress(0.5,detail="Erasing Temporal Files")
        file.remove(paste0('./',input$file1$name))
        incProgress(1,detail = 'Done')
      })
      removeModal()
      showModal(modalDialog('File Annotated',
                            p('File Available at: ',paste0(cfg$PATHS$resultspath,'vDBannotated/',gsub('\\b.xlsx|\\b.annotated','',input$file1$name),'.annotated.xlsx')
                            ),
                            downloadButton('downloadannotatedVDB','Download Result')
                            )


                )
    }
    })
  observeEvent(input$addbtn, {
    if((input$url=="")&(is.null(input$file1))){
      removeModal()
      modalDialog('No input data to download.')
    }
    else if((input$url=="")&!(is.null(input$file1))){
      system2('MODApy',args = paste('addPatient',input$file1$datapath),wait = FALSE,stdout = FALSE,stderr = FALSE)
      removeModal()
      showModal(modalDialog('Downloading'))
    }
    else if(!(input$url=="")&(is.null(input$file1))){
      system2('MODApy',args = paste('addPatient',input$url),wait = FALSE,stdout = FALSE,stderr = FALSE)
      removeModal()
      showModal(modalDialog('Downloading'))
    }
    else if(!(input$url=="")&!(is.null(input$file1))){
      system2('MODApy',args = paste('addPatient',input$file1$datapath),wait = FALSE,stdout = FALSE,stderr = FALSE)
      removeModal()
      showModal(modalDialog('Downloading'))
    }

  })
  observeEvent(input$buttonrun, {
    shinyjs::disable('buttonrun')
    file.create(logfile)
    rv$started<-TRUE
    getcommand(input)
  })
  ##Runs Pipeline for Fast Q
  observeEvent(input$runPipelinefq, {
    #shinyjs::disable('buttonrun')
    patpath = gsub('_1','',input$fqfile)
    if(length(system('ps aux | grep "MODApy pipeline"',intern=TRUE))>2){
      showModal(pipelineRunningModal())
    }
    else if(file.exists(paste(cfg$PATHS$patientpath, patpath, "/",patpath,"_1", ".fastq", sep=""))){
      file.create(logfile)
      cmd =  paste("pipeline -Pipeline",
                   "BestPractices-Trim.json", "-FQ",
                   paste(patpath, "/",patpath,"_1", ".fastq", sep=""), "-FQ",
                   paste(patpath, "/",patpath,"_2", ".fastq", sep=""))
      print(cmd)
      system2("MODApy", cmd ,wait=FALSE,stdout = FALSE,stderr = FALSE)

    }
    else if(file.exists(paste(cfg$PATHS$patientpath, patpath, "/",patpath,"_1", ".fastq.gz", sep=""))){
      file.create(logfile)
      cmd = paste("pipeline -Pipeline",
                  "BestPractices-Trim.json", "-FQ",
                  paste(patpath, "/",patpath,"_1", ".fastq.gz", sep=""), "-FQ",
                  paste(patpath, "/",patpath,"_2", ".fastq.gz", sep=""))
      system2("MODApy", cmd ,wait=FALSE,stdout = FALSE,stderr = FALSE)

    }
    else if(file.exists(paste(cfg$PATHS$patientpath, patpath, "/",patpath,"_1", ".fq.gz", sep=""))){
      file.create(logfile)
      cmd = paste("pipeline -Pipeline",
                  "BestPractices-Trim.json", "-FQ",
                  paste(patpath, "/",patpath,"_1", ".fq.gz", sep=""), "-FQ",
                  paste(patpath, "/",patpath,"_2", ".fq.gz", sep=""))
      system2("MODApy", cmd ,wait=FALSE,stdout = FALSE,stderr = FALSE)

    }
    else if(file.exists(paste(cfg$PATHS$patientpath, patpath, "/",patpath,"_1", ".fq", sep=""))){
      file.create(logfile)
      cmd = paste("pipeline -Pipeline",
                  "BestPractices-Trim.json", "-FQ",
                  paste(patpath, "/",patpath,"_1", ".fq", sep=""), "-FQ",
                  paste(patpath, "/",patpath,"_2", ".fq", sep=""))
      system2("MODApy", cmd ,wait=FALSE,stdout = FALSE,stderr = FALSE)
    }
  })
  ##Runs Pipeline for Bam
  observeEvent(input$runPipelinebam, {
    #shinyjs::disable('buttonrun')
    patpath=paste0(cfg$PATHS$patientpath,input$bamfile,'/',input$bamfile,'.bam')
    if(length(system('ps aux | grep "MODApy pipeline"',intern=TRUE))>2){
      showModal(pipelineRunningModal())
    }
    else if(file.exists(patpath)){
      cmd =  paste("pipeline -Pipeline",
                   "BestPracticesBam-Annotated.json", "-FQ",
                   paste0(input$bamfile,'/',input$bamfile,'.bam'))
      file.create(logfile)
      system2("MODApy", cmd ,wait=FALSE,stdout = FALSE,stderr = FALSE)
    }})
    ##Runs Pipeline for VCF
  observeEvent(input$runPipelinevcf, {
    #shinyjs::disable('buttonrun')
    patpath=paste0(cfg$PATHS$patientpath,input$vcffile,'/',input$vcffile,'.final.vcf')
    if(length(system('ps aux | grep "MODApy pipeline"',intern=TRUE))>2){
      showModal(pipelineRunningModal())
    }
    else if(file.exists(patpath)){
      file.create(logfile)
      cmd =  paste("pipeline -Pipeline",
                    "BestPractices-Annotation.json", "-FQ",
                    paste0(input$vcffile,'/',input$vcffile,'.final.vcf'))
      system2("MODApy", cmd ,wait=FALSE,stdout = FALSE,stderr = FALSE)
  }})
  #Download Result Annotated
  output$downloadannotatedVDB <- downloadHandler(
    filename <- paste0(cfg$PATHS$resultspath,'vDBannotated/',gsub('\\b.xlsx|\\b.annotated','',input$file1$name),'.annotated.xlsx'),
    content<-function(downfile){
      filepath <- paste0(cfg$PATHS$resultspath,'vDBannotated/',gsub('\\b.xlsx|\\b.annotated','',input$file1$name),'.annotated.xlsx')
      print(filepath)
      file.copy(filepath,downfile)
      file.size(filepath)
    },
    contentType = "application/xlsx"
  )
  #Download Result
  output$downloadData <- downloadHandler(
    filename<-function(){
      downpath$file
    },
    content<-function(file){
      file.copy(downpath$path,file)
    },
    contentType = "application/xlsx"
  )
  observeEvent(input$buttonlastcmd, {
    rv$textstream = ""
    rv$started<-TRUE
  })
  observeEvent(
    input$Patient1D,{
      current <- isolate(input$Patient2D)
      updateSelectInput(session,'Patient2D',selected = current, choices = patientsvcf[!(patientsvcf %in% input$Patient1D)])
    })
  observeEvent(
    input$Patient2D,{
      current <- isolate(input$Patient1D)
      updateSelectInput(session,'Patient1D',selected = current, choices = patientsvcf[!(patientsvcf %in% input$Patient2D)])
    }
  )
  observeEvent(
    input$Patient1T,{
      current2 <- isolate(input$Patient2T)
      current3 <- isolate(input$Patient3T)
      updateSelectInput(session,'Patient2T',selected = current2,choices = patientsvcf[!(patientsvcf %in% c(input$Patient1T,input$Patient3T))])
      updateSelectInput(session,'Patient3T',selected = current3,choices = patientsvcf[!(patientsvcf %in% c(input$Patient1T,input$Patient2T))])
    })
  observeEvent(
    input$Patient2T,{
      current1 <- isolate(input$Patient1T)
      current3 <- isolate(input$Patient3T)
      updateSelectInput(session,'Patient1T',selected = current1,choices = patientsvcf[!(patientsvcf %in% c(input$Patient2T,input$Patient3T))])
      updateSelectInput(session,'Patient3T',selected = current3,choices = patientsvcf[!(patientsvcf %in% c(input$Patient1T,input$Patient2T))])
    })
  observeEvent(
    input$Patient3T,{
      current1 <- isolate(input$Patient1T)
      current2 <- isolate(input$Patient2T)
      updateSelectInput(session,'Patient1T',selected = current1,choices = patientsvcf[!(patientsvcf %in% c(input$Patient2T,input$Patient3T))])
      updateSelectInput(session,'Patient2T',selected = current2,choices = patientsvcf[!(patientsvcf %in% c(input$Patient1T,input$Patient3T))])
    })
  observe({
    rv$timer()
    if(isolate(rv$started))rv$textstream <- paste(readLines(logfile),collapse="<br/>")
    if(grepl('Complete',rv$textstream,ignore.case = TRUE)){
      rv$started<-FALSE
      dwnpath <- strsplit(rv$textstream,'File available at:')[[1]][2]
      dwnname <- tail(unlist(strsplit(strsplit(rv$textstream,'File available at:')[[1]][2],'/')),n=1)
      downpath$path <- dwnpath
      downpath$file <- dwnname
      rv$textstream <- gsub('Complete','Finished',rv$textstream,fixed = TRUE)
      shinyjs::enable('buttonrun')
      shinyjs::enable('downloadData')
    }
    else if(grepl('Failed',rv$textstream,ignore.case = TRUE)){
      rv$started<-FALSE
      rv$textstream <- gsub('Failed','had an error.',rv$textstream,fixed = TRUE)
      shinyjs::enable('buttonrun')
    }
  })

  observe({
    rv2$timer()
    Sys.sleep(0.2)
    if(isolate(rv2$started))rv$textstream2 <- gsub('\\{|\\}|\\"|\\,','',paste(readLines(dlog,warn=FALSE),collapse="<br/>"))
  })
  observe({
    rv3$timer()
    Sys.sleep(0.2)
    if(isolate(rv3$started))rv$textstream3 <- paste(readLines(pipelog),collapse="<br/>")
  })
  output$downout <- renderUI({
    HTML(rv$textstream2)
  })
  output$consoleout <- renderUI({
    HTML(rv$textstream)
  })
  output$pipeout <- renderUI({
    HTML(rv$textstream3)
  })
  output$dbout <- renderUI({
    HTML(rv$textstream)
  })
  output$vennDuos <- renderUI({
    selectInput("vennplaceD",NULL, choices = list(input$Patient1D, input$Patient2D, paste(input$Patient1D, input$Patient2D,sep=":"), "All"), selected='All')
  })
  output$vennTrios <- renderUI({
    selectInput("vennplaceT",NULL, choices = list(input$Patient1T, input$Patient2T, input$Patient3T, paste(input$Patient1T, input$Patient2T,sep=":"), paste(input$Patient1T, input$Patient3T,sep=":"), paste(input$Patient2T, input$Patient3T,sep=":"), paste(input$Patient1T, input$Patient2T, input$Patient3T,sep=":"), "All"), selected='All')
  })
  observeEvent(input$cfgButton,{
    showModal(configModal())
  })
  observeEvent(input$savecfg,{
    print(typeof(cfg$GENERAL))
  })

  session$onSessionEnded(function(){stopApp()})
}

shinyApp(ui = ui, server = server)
