library(shiny)
library(ConfigParser)
library(reticulate)
library(DT)
library(openxlsx)
# python config -------------------------------------------------------------------
cfgpath = '/DiscoDatos/Development/modapy/MODApy/config.ini'
logfile = "/DiscoDatos/Development/modapy/MODApy/logs/currentrun.log"
dlog = "/DiscoDatos/Development/modapy/MODApy/logs/downloads.log"
cfg = read.ini(cfgpath)
use_virtualenv("/DiscoDatos/Development/modapy/venv/")
use_python("/DiscoDatos/Development/modapy/venv/python3")
py_config()


# combo box options -------------------------------------------------------------------
patientsvcf <- gsub('\\..*','',basename(list.files(cfg$PATHS$patientpath,pattern="\\.final.vcf",recursive = TRUE)))

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
          Panel={
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
panels =
  ui <- navbarPage("MODApy",
                   tabPanel("Analisis",
                        sidebarLayout(
                          sidebarPanel(
                            width = 6,
                            tabsetPanel(id="tabset", selected = "Single",
                                        # tabPanel("Pipeline",
                                        #   selectInput(inputId = "Pipeline", label = "Pipelines", choices = list.files(path=cfg$PATHS$pipelinespath)),
                                        #   selectInput(inputId = "PatientPipe", label = "Patient", choices = list.dirs(
                                        #     path = cfg$PATHS$patientpath ,full.names = FALSE,recursive = FALSE))
                                        # ),
                                        tabPanel("Single",
                                                 selectInput(inputId = "PatientPanel", label = "Patient", choices = patientsvcf),
                                                 selectInput(inputId = "Panel", label = "Panel", choices = gsub('.xlsx','',list.files(path=cfg$PATHS$panelspath)))
                                        ),
                                        tabPanel("Duos",
                                                 selectInput(inputId = "Patient1D", label = "Patient 1", choices = patientsvcf),
                                                 selectInput(inputId = "Patient2D", label = "Patient 2", choices = patientsvcf),
                                                 selectInput(inputId = "PanelD", label = "Panel (optional)", choices = gsub('.xlsx',"",c('NONE',list.files(path=cfg$PATHS$panelspath))), selected='NONE'),
                                                 uiOutput("vennDuos")
                                        ),
                                        tabPanel("Trios",
                                                 selectInput(inputId = "Patient1T", label = "Patient 1", choices = patientsvcf),
                                                 selectInput(inputId = "Patient2T", label = "Patient 2", choices = patientsvcf),
                                                 selectInput(inputId = "Patient3T", label = "Patient 3", choices = patientsvcf),
                                                 selectInput(inputId = "PanelT", label = "Panel (optional)", choices = gsub('.xlsx',"",c('NONE',list.files(path=cfg$PATHS$panelspath))), selected='NONE'),
                                                 uiOutput("vennTrios")
                                        )
                            ),
                            actionButton("buttonrun","Run"),
                            actionButton("buttonlastcmd","Get Last Command Status")
                          ),
                          mainPanel(#style="background-color:blue",
                            width = 6,
                            h1('Command Output',style = "font-family: 'Courgette', cursive;
                               font-weight: 500; line-height: 1.1; 
                               color: #4d3a7d;"),
                            htmlOutput("consoleout")
                            )
                          )),
               tabPanel('Add New Patient',
                        p('This will add a new patient, downlading the data from a url or a xlsx/xls file.
                          Please either write the url or upload a file and press Submit. If both values are filled, will only download urls from files. 
                          It will generate the folder under the Patient folder and download the .tar file'),
                        textInput('url',NULL,value="",placeholder='Enter URL'),
                        fileInput('file1','Choose File to Upload',accept=c('.xls','.xlsx')),
                        actionButton("addbtn","Add Patient"),
                        actionButton("dstatbtn","Get Downloads Status"),
                        h1('Command Output',style = "font-family: 'Courgette', cursive;
                           font-weight: 500; line-height: 1.1; 
                           color: #4d3a7d;"),
                        htmlOutput("downout")
                        ),
               tabPanel('VariantsDB',
                        actionButton("buildDBbtn","Build DataBase"),
                        actionButton("openDBbtn","Open Database"),
                        htmlOutput("dbout"),
                        DT::dataTableOutput("mytable")
                        )
    )
# server -------------------------------------------------------------------
server <- function(input,output, session){
  rv <- reactiveValues(textstream = c(""), timer = reactiveTimer(1000),started=FALSE)
  rv2 <- reactiveValues(textstream2 = c(""), timer = reactiveTimer(1000),started=FALSE)
  rv3 <- reactiveValues(textstream2 = c(""), timer = reactiveTimer(1000),started=FALSE)

  observeEvent(input$buildDBbtn, {
    rv$textstream = ""
    rv$started<-TRUE
    system2('MODApy',args = 'variantsDB -buildDB',wait = FALSE,stdout = FALSE,stderr = FALSE)
  })
  observeEvent(input$openDBbtn, {
    if(file.exists(cfg$PATHS$dbpath)){
      rv$textstream = ""
      rv$started<-FALSE
      #df1 <-read.xlsx(cfg$PATHS$dbpath, sheet=1,skipEmptyRows=FALSE)
      df1 <- read.csv(cfg$PATHS$dbpath,check.names = FALSE)
      output$mytable = DT::renderDataTable({df1})
    }
    else{
      rv$textstream = "Variants file not found. Try to build variantsDB first."
    }

  })
  observeEvent(input$addbtn, {
    if((input$url=="")&(is.null(input$file1))){
      rv2$started<-FALSE
      rv3$started<-FALSE
      rv$textstream2 <- "No input selected to download."
    }
    else if((input$url=="")&!(is.null(input$file1))){
      rv2$started<-FALSE
      rv3$started<-TRUE
      system2('MODApy',args = paste('addPatient',input$file1$datapath),wait = FALSE,stdout = FALSE,stderr = FALSE)
    }
    else if(!(input$url=="")&(is.null(input$file1))){
      rv2$started<-FALSE
      rv3$started<-TRUE
      system2('MODApy',args = paste('addPatient',input$url),wait = FALSE,stdout = FALSE,stderr = FALSE)
    }
    else if(!(input$url=="")&!(is.null(input$file1))){
      rv2$started<-FALSE
      rv3$started<-TRUE
      system2('MODApy',args = paste('addPatient',input$file1$datapath),wait = FALSE,stdout = FALSE,stderr = FALSE)}
  })
  observeEvent(input$buttonrun, {
    rv$textstream = ""
    rv$started<-TRUE
    getcommand(input)
  })
  observeEvent(input$buttonlastcmd, {
    rv$textstream = ""
    rv$started<-TRUE
  })
  observeEvent(input$dstatbtn, {
    rv3$started<-FALSE
    rv2$started<-TRUE
  })
  observe({
    rv$timer()
    Sys.sleep(0.2)
    if(isolate(rv$started))rv$textstream <- paste(readLines(logfile),collapse="<br/>")
  })
  observe({
    rv2$timer()
    Sys.sleep(0.2)
    if(isolate(rv2$started))rv$textstream2 <- gsub('\\{|\\}|\\"|\\,','',paste(readLines(dlog,warn=FALSE),collapse="<br/>"))
  })
  observe({
    rv3$timer()
    Sys.sleep(0.2)
    if(isolate(rv3$started))rv$textstream2 <- paste(readLines(logfile),collapse="<br/>")
  })
  output$downout <- renderUI({
    HTML(rv$textstream2)
  })
  output$consoleout <- renderUI({
    HTML(rv$textstream)
  })
  output$dbout <- renderUI({
    HTML(rv$textstream)
  })
  output$vennDuos <- renderUI({
    selectInput("vennplaceD","Venn Place:", choices = list(input$Patient1D, input$Patient2D, paste(input$Patient1D, input$Patient2D,sep=":"), "All"), selected='All')
  })
  output$vennTrios <- renderUI({
    selectInput("vennplaceT","Venn Place:", choices = list(input$Patient1T, input$Patient2T, input$Patient3T, paste(input$Patient1T, input$Patient2T,sep=":"), paste(input$Patient1T, input$Patient3T,sep=":"), paste(input$Patient2T, input$Patient3T,sep=":"), paste(input$Patient1T, input$Patient2T, input$Patient3T,sep=":"), "All"), selected='All')
  })
  session$onSessionEnded(function(){stopApp()})
}

shinyApp(ui = ui, server = server)