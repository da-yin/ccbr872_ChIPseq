library(shiny)
library(plotly)
library(ggplot2)
library(dplyr)
library(plyr)

#options(shiny.port=36569)
#data_path = '/data/CCBR_Pipeliner/db/PipeDB/db/GTRD/'
data_path = '~/Desktop/active_projects/ccbr872_ChIPseq/'


file_suffix = '.intersummit_distance.txt'

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Intersummit Distance"),
  # 
  # tags$head(
  #     tags$script(src = "https://cdnjs.cloudflare.com/ajax/libs/jQuery.print/1.6.0/jQuery.print.min.js")
  # ),
  
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      textInput("TF1", "TF1", "Hnf1b"),
      textInput("TF2", "TF2", "Chd1"),
      textInput("TF3", "TF3", "Mafa"),
      textInput("TF4", "TF4", "Mllt3"),
      br(),
      actionButton("boxplotButton", "Generate Boxplots"),
      br(),
      br(),
      actionButton("densityButton", "Generate Density Plot"),
      br(),
      br(),
      actionButton("ksButton", "Kolmogorov-Smirnov test"),
      br(),
      br(),
      #actionButton("statsButton", "generate stats"),
      br(),
      br(),
      # actionButton("print", "Print", onclick = "$('#textarea').print();")
    ),
    
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("boxplots"),
      plotOutput("densityplots"),
      verbatimTextOutput("stats1"),
      verbatimTextOutput("ksTest1"),
      verbatimTextOutput("stats2"),
      verbatimTextOutput("ksTest2"),
      verbatimTextOutput("stats3"),
      verbatimTextOutput("ksTest3")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  
  get_genes <- eventReactive(input$boxplotButton, {
    # Validation. All 4 genes must be typed
    # It also validates if the files exist and they're readable
    validate(
      need(input$TF1 != "", "Please type a valid gene 1"),
      need(input$TF2 != "", "Please type a valid gene 2"),
      need(input$TF3 != "", "Please type a valid gene 3"),
      need(input$TF4 != "", "Please type a valid gene 4")
      # need(try(read.table(paste0(data_path, input$TF1,'_', input$TF2, file_suffix), sep='\t', header=FALSE)), 
      #      paste0("There is no file for the pair ", input$TF1, ' and ', input$TF2)),
      # need(try(read.table(paste0(data_path, input$TF3,'_', input$TF4, file_suffix), sep='\t', header=FALSE)), 
      #      paste0("There is no file for the pair ", input$TF3, ' and ', input$TF4))
    )
    
    c(input$TF1, input$TF2, input$TF3, input$TF4)
  })
  
  
  get_genes_density <- eventReactive(input$densityButton, {
    # Validation. All 4 genes must be typed
    # It also validates if the files exist and they're readable
    validate(
      need(input$TF1 != "", "Please type a valid gene 1"),
      need(input$TF2 != "", "Please type a valid gene 2"),
      need(input$TF3 != "", "Please type a valid gene 3"),
      need(input$TF4 != "", "Please type a valid gene 4")
      # need(try(read.table(paste0(data_path, input$TF1,'_', input$TF2, file_suffix), sep='\t', header=FALSE)), 
      #      paste0("There is no file for the pair ", input$TF1, ' and ', input$TF2)),
      # need(try(read.table(paste0(data_path, input$TF3,'_', input$TF4, file_suffix), sep='\t', header=FALSE)), 
      #      paste0("There is no file for the pair ", input$TF1, ' and ', input$TF2))
    )
    
    c(input$TF1, input$TF2, input$TF3, input$TF4)
  })
  
  get_genes_ks <- eventReactive(input$ksButton, {
    # Validation. All 4 genes must be typed
    # It also validates if the files exist and they're readable
    validate(
      need(input$TF1 != "", "Please type a valid gene 1"),
      need(input$TF2 != "", "Please type a valid gene 2"),
      need(input$TF3 != "", "Please type a valid gene 3"),
      need(input$TF4 != "", "Please type a valid gene 4")
      # need(try(read.table(paste0(data_path, input$TF1,'_', input$TF2, file_suffix), sep='\t', header=FALSE)), 
      #      paste0("There is no file for the pair ", input$TF1, ' and ', input$TF2)),
      # need(try(read.table(paste0(data_path, input$TF3,'_', input$TF4, file_suffix), sep='\t', header=FALSE)), 
      #      paste0("There is no file for the pair ", input$TF1, ' and ', input$TF2))
    )
    
    c(input$TF1, input$TF2, input$TF3, input$TF4)
  })
  
  output$boxplots <- renderPlot({
    genes = get_genes()
    a = paste0(data_path, genes[1],'_', genes[2], file_suffix)
    b = paste0(data_path, genes[2],'_', genes[1], file_suffix)
    c = paste0(data_path, genes[3],'_', genes[4], file_suffix)
    d = paste0(data_path, genes[4],'_', genes[3], file_suffix)
    
    if (file.exists(a)){
      file1 = a
    }
    if (file.exists(b)){
      file1 = b
    }
    if (file.exists(c)){
      file2 = c
    }
    if (file.exists(d)){
      file2 = d
    }
    
    # file1 = paste0(data_path, genes[1],'_', genes[2], file_suffix)
    # file2 = paste0(data_path, genes[3],'_', genes[4], file_suffix)
    data1 <- read.table(file1, sep='\t', header=FALSE)
    data2 <- read.table(file2, sep='\t', header=FALSE)
    
    data1$pair = paste0(data1$V1,"_",data1$V2)
    data2$pair = paste0(data2$V1,"_",data2$V2)
    
    # control for genes1 and genes2 
    my_files1 = list.files(pattern = genes[1])
    my_files2 = list.files(pattern = genes[2])
    
    my_distance1 <- lapply(my_files1, function(x) {
      tryCatch(read.table(x, header = FALSE, sep = "\t"), error=function(e) NULL)
    })
    
    my_distance2 <- lapply(my_files2, function(x) {
      tryCatch(read.table(x, header = FALSE, sep = "\t"), error=function(e) NULL)
    })
    
    data3.1 = do.call(rbind.fill, my_distance1)
    data3.2 = do.call(rbind.fill, my_distance2)
    data3 = rbind(data3.1, data3.2)
    data3$pair = paste0(genes[2],'_', genes[1],"_control")
    
    # control for genes3 and genes4
    my_files3 = list.files(pattern = genes[3])
    my_files4 = list.files(pattern = genes[4])
    
    my_distance3 <- lapply(my_files3, function(x) {
      tryCatch(read.table(x, header = FALSE, sep = "\t"), error=function(e) NULL)
    })
    
    my_distance4 <- lapply(my_files4, function(x) {
      tryCatch(read.table(x, header = FALSE, sep = "\t"), error=function(e) NULL)
    })
    
    data4.1 = do.call(rbind.fill, my_distance3)
    data4.2 = do.call(rbind.fill, my_distance4)
    data4 = rbind(data4.1, data4.2)
    data4$pair = paste0(genes[4],'_', genes[3],"_control")
    
    
    df_both = rbind(data1, data2, data3, data4)
    colnames(df_both)[3]="distance"
    df_both = df_both[,c(3,4)]
    df_both$distance=abs(df_both$distance)
    df_both$distance[mapply(is.infinite, df_both$distance)] = 0
    theme_set(theme_classic())
    qplot(pair, distance, data = df_both, 
          geom=c("boxplot"), fill = pair,log="y")
  })
  
  
  output$densityplots <- renderPlot({
    genes = get_genes_density()
    
    a = paste0(data_path, genes[1],'_', genes[2], file_suffix)
    b = paste0(data_path, genes[2],'_', genes[1], file_suffix)
    c = paste0(data_path, genes[3],'_', genes[4], file_suffix)
    d = paste0(data_path, genes[4],'_', genes[3], file_suffix)
    
    if (file.exists(a)){
      file1 = a
    }
    if (file.exists(b)){
      file1 = b
    }
    if (file.exists(c)){
      file2 = c
    }
    if (file.exists(d)){
      file2 = d
    }
    
    
    # file1 = paste0(data_path, genes[1],'_', genes[2], file_suffix)
    # file2 = paste0(data_path, genes[3],'_', genes[4], file_suffix)
    data1 <- read.table(file1, sep='\t', header=FALSE)
    
    data2 <- read.table(file2, sep='\t', header=FALSE)
    
    data1$pair = paste0(data1$V1,"_",data1$V2)
    data2$pair = paste0(data2$V1,"_",data2$V2)
    
    my_files1 = list.files(pattern = genes[1])
    my_files2 = list.files(pattern = genes[2])
    my_files3 = list.files(pattern = genes[3])
    my_files4 = list.files(pattern = genes[4])
    
    my_distance1 <- lapply(my_files1, function(x) {
      tryCatch(read.table(x, header = FALSE, sep = "\t"), error=function(e) NULL)
    })
    
    my_distance2 <- lapply(my_files2, function(x) {
      tryCatch(read.table(x, header = FALSE, sep = "\t"), error=function(e) NULL)
    })
    
    my_distance3 <- lapply(my_files3, function(x) {
      tryCatch(read.table(x, header = FALSE, sep = "\t"), error=function(e) NULL)
    })
    
    my_distance4 <- lapply(my_files4, function(x) {
      tryCatch(read.table(x, header = FALSE, sep = "\t"), error=function(e) NULL)
    })
    
    
    data3.1 = do.call(rbind.fill, my_distance1)
    data3.2 = do.call(rbind.fill, my_distance2)
    data3 = rbind(data3.1, data3.2)
    data3$pair = paste0(genes[2],'_', genes[1],"_control")
    data4.1 = do.call(rbind.fill, my_distance3)
    data4.2 = do.call(rbind.fill, my_distance4)
    data4 = rbind(data4.1, data4.2)
    data4$pair = paste0(genes[4],'_', genes[3],"_control")
    
    df_both = rbind(data1, data2, data3, data4)
    colnames(df_both)[3]="distance"
    df_both = df_both[,c(3,4)]
    df_both$distance=abs(df_both$distance)
    df_both$distance[mapply(is.infinite, df_both$distance)] = 0
    
    options(scipen=10000)
    theme_set(theme_classic())
    qplot(distance, data = df_both, geom = "density", color = pair,y = ..scaled..,log="x")
  })
  
  output$stats1 <- renderPrint({
    genes = get_genes_ks()
    a = paste0(data_path, genes[1],'_', genes[2], file_suffix)
    b = paste0(data_path, genes[2],'_', genes[1], file_suffix)
    if (file.exists(a)){
      file1 = a
    }
    if (file.exists(b)){
      file1 = b
    }
    c = paste0(data_path, genes[3],'_', genes[4], file_suffix)
    d = paste0(data_path, genes[4],'_', genes[3], file_suffix)
    if (file.exists(c)){
      file2 = c
    }
    if (file.exists(d)){
      file2 = d
    }
    
    data1 <- read.table(file1, sep='\t', header=FALSE)
    abs_data1 = abs(data1[3])
    
    
    data2 <- read.table(file2, sep='\t', header=FALSE)
    abs_data2 = abs(data2[3])
    
    paste0("Median intersummit distance ", genes[1],'_', genes[2], " : ",median(unlist(abs_data1)),
           "; ", genes[3],'_', genes[4], " : ",median(unlist(abs_data2)))
    
  })
  
  output$ksTest3 <- renderPrint({
    genes = get_genes_ks()
    c = paste0(data_path, genes[3],'_', genes[4], file_suffix)
    d = paste0(data_path, genes[4],'_', genes[3], file_suffix)
    if (file.exists(c)){
      file2 = c
    }
    if (file.exists(d)){
      file2 = d
    }    
    data2 <- read.table(file2, sep='\t', header=FALSE)
    TF3_TF4 = unlist(abs(data2[3]))
    
    # control for genes3 and genes4
    my_files3 = list.files(pattern = genes[3])
    my_files4 = list.files(pattern = genes[4])
    
    my_distance3 <- lapply(my_files3, function(x) {
      tryCatch(read.table(x, header = FALSE, sep = "\t"), error=function(e) NULL)
    })
    
    my_distance4 <- lapply(my_files4, function(x) {
      tryCatch(read.table(x, header = FALSE, sep = "\t"), error=function(e) NULL)
    })
    
    data4.1 = do.call(rbind.fill, my_distance3)
    data4.2 = do.call(rbind.fill, my_distance4)
    data4 = rbind(data4.1, data4.2)
    data4$pair = paste0(genes[4],'_', genes[3],"_control")
    abs_data4 = abs(data4[3])
    
    TF3_TF4_control = unlist(abs(data4[3]))
    
    ks = ks.test(TF3_TF4,TF3_TF4_control)
    
    ks
  })
  
  output$ksTest2 <- renderPrint({
    genes = get_genes_ks()
    
    a = paste0(data_path, genes[1],'_', genes[2], file_suffix)
    b = paste0(data_path, genes[2],'_', genes[1], file_suffix)
    if (file.exists(a)){
      file1 = a
    }
    if (file.exists(b)){
      file1 = b
    }
    data1 <- read.table(file1, sep='\t', header=FALSE)
    TF1_TF2 = unlist(abs(data1[3]))
    
    
    my_files1 = list.files(pattern = genes[1])
    my_files2 = list.files(pattern = genes[2])
    
    my_distance1 <- lapply(my_files1, function(x) {
      tryCatch(read.table(x, header = FALSE, sep = "\t"), error=function(e) NULL)
    })
    
    my_distance2 <- lapply(my_files2, function(x) {
      tryCatch(read.table(x, header = FALSE, sep = "\t"), error=function(e) NULL)
    })
    
    data3.1 = do.call(rbind.fill, my_distance1)
    data3.2 = do.call(rbind.fill, my_distance2)
    data3 = rbind(data3.1, data3.2)
    data3$pair = paste0(genes[2],'_', genes[1],"_control")
    TF1_TF2_control = unlist(abs(data3[3]))
    
    
    # ks = ks.test(unlist(data1[3]), unlist(data2[3]))
    ks = ks.test(TF1_TF2,TF1_TF2_control)
    
    ks
  })
  
  
  output$stats3 <- renderPrint({
    genes = get_genes_ks()
    c = paste0(data_path, genes[3],'_', genes[4], file_suffix)
    d = paste0(data_path, genes[4],'_', genes[3], file_suffix)
    if (file.exists(c)){
      file2 = c
    }
    if (file.exists(d)){
      file2 = d
    }
    data2 <- read.table(file2, sep='\t', header=FALSE)
    abs_data2 = abs(data2[3])
    
    # control for genes3 and genes4
    my_files3 = list.files(pattern = genes[3])
    my_files4 = list.files(pattern = genes[4])
    
    my_distance3 <- lapply(my_files3, function(x) {
      tryCatch(read.table(x, header = FALSE, sep = "\t"), error=function(e) NULL)
    })
    
    my_distance4 <- lapply(my_files4, function(x) {
      tryCatch(read.table(x, header = FALSE, sep = "\t"), error=function(e) NULL)
    })
    
    data4.1 = do.call(rbind.fill, my_distance3)
    data4.2 = do.call(rbind.fill, my_distance4)
    data4 = rbind(data4.1, data4.2)
    data4$pair = paste0(genes[4],'_', genes[3],"_control")
    abs_data4 = abs(data4[3])
    
    
    paste0("Median intersummit distance ", genes[3],'_', genes[4], " : ",median(unlist(abs_data2)),
           "; ", genes[4],'_', genes[3],"_control", " : ",median(unlist(abs_data4)))
    
    
  })
  
  output$stats2 <- renderPrint({
    genes = get_genes_ks()
    a = paste0(data_path, genes[1],'_', genes[2], file_suffix)
    b = paste0(data_path, genes[2],'_', genes[1], file_suffix)
    if (file.exists(a)){
      file1 = a
    }
    if (file.exists(b)){
      file1 = b
    }
    data1 <- read.table(file1, sep='\t', header=FALSE)
    abs_data1 = abs(data1[3])
    
    my_files1 = list.files(pattern = genes[1])
    my_files2 = list.files(pattern = genes[2])
    
    my_distance1 <- lapply(my_files1, function(x) {
      tryCatch(read.table(x, header = FALSE, sep = "\t"), error=function(e) NULL)
    })
    
    my_distance2 <- lapply(my_files2, function(x) {
      tryCatch(read.table(x, header = FALSE, sep = "\t"), error=function(e) NULL)
    })
    
    data3.1 = do.call(rbind.fill, my_distance1)
    data3.2 = do.call(rbind.fill, my_distance2)
    data3 = rbind(data3.1, data3.2)
    data3$pair = paste0(genes[2],'_', genes[1],"_control")
    abs_data3 = abs(data3[3])
    
    paste0("Median intersummit distance ", genes[1],'_', genes[2], " : ",median(unlist(abs_data1)),
           "; ", genes[2],'_', genes[1],"_control", " : ",median(unlist(abs_data3)))
    
    
  })
  
  
  
  output$ksTest1 <- renderPrint({
    genes = get_genes_ks()
    a = paste0(data_path, genes[1],'_', genes[2], file_suffix)
    b = paste0(data_path, genes[2],'_', genes[1], file_suffix)
    if (file.exists(a)){
      file1 = a
    }
    if (file.exists(b)){
      file1 = b
    }
    c = paste0(data_path, genes[3],'_', genes[4], file_suffix)
    d = paste0(data_path, genes[4],'_', genes[3], file_suffix)
    if (file.exists(c)){
      file2 = c
    }
    if (file.exists(d)){
      file2 = d
    }
    data1 <- read.table(file1, sep='\t', header=FALSE)
    TF1_TF2 = unlist(abs(data1[3]))
    
    data2 <- read.table(file2, sep='\t', header=FALSE)
    TF3_TF4 = unlist(abs(data2[3]))
    
    # Replace -Inf caused by log-zero values
    TF1_TF2[mapply(is.infinite, TF1_TF2)] = 0
    TF3_TF4[mapply(is.infinite, TF3_TF4)] = 0
    
    
    # ks = ks.test(unlist(data1[3]), unlist(data2[3]))
    ks = ks.test(TF1_TF2, TF3_TF4)
    
    ks
  })
}




# Run the application 
shinyApp(ui = ui, server = server)
