library(shiny)
library(plotly)
library(ggplot2)
library(dplyr)

data_path = '~/Desktop/active_projects/ccbr872_ChIPseq/rawdata/'
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
            textInput("gene1", "gene1", "Hnf1b"),
            textInput("gene2", "gene2", "Chd1"),
            textInput("gene3", "gene3", "Mafa"),
            textInput("gene4", "gene4", "Mllt3"),
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
            verbatimTextOutput("stats2"),
            verbatimTextOutput("ksTest")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    
    get_genes <- eventReactive(input$boxplotButton, {
        # Validation. All 4 genes must be typed
        # It also validates if the files exist and they're readable
        validate(
            need(input$gene1 != "", "Please type a valid gene 1"),
            need(input$gene2 != "", "Please type a valid gene 2"),
            need(input$gene3 != "", "Please type a valid gene 3"),
            need(input$gene4 != "", "Please type a valid gene 4"),
            need(try(read.table(paste0(data_path, input$gene1,'_', input$gene2, file_suffix), sep='\t', header=FALSE)), 
                 paste0("There is no file for the pair ", input$gene1, ' and ', input$gene2)),
            need(try(read.table(paste0(data_path, input$gene3,'_', input$gene4, file_suffix), sep='\t', header=FALSE)), 
                 paste0("There is no file for the pair ", input$gene3, ' and ', input$gene4))
        )
        
        c(input$gene1, input$gene2, input$gene3, input$gene4)
    })
    
    
    get_genes_density <- eventReactive(input$densityButton, {
        # Validation. All 4 genes must be typed
        # It also validates if the files exist and they're readable
        validate(
            need(input$gene1 != "", "Please type a valid gene 1"),
            need(input$gene2 != "", "Please type a valid gene 2"),
            need(input$gene3 != "", "Please type a valid gene 3"),
            need(input$gene4 != "", "Please type a valid gene 4"),
            need(try(read.table(paste0(data_path, input$gene1,'_', input$gene2, file_suffix), sep='\t', header=FALSE)), 
                 paste0("There is no file for the pair ", input$gene1, ' and ', input$gene2)),
            need(try(read.table(paste0(data_path, input$gene3,'_', input$gene4, file_suffix), sep='\t', header=FALSE)), 
                 paste0("There is no file for the pair ", input$gene1, ' and ', input$gene2))
        )
        
        c(input$gene1, input$gene2, input$gene3, input$gene4)
    })
    
    get_genes_ks <- eventReactive(input$ksButton, {
        # Validation. All 4 genes must be typed
        # It also validates if the files exist and they're readable
        validate(
            need(input$gene1 != "", "Please type a valid gene 1"),
            need(input$gene2 != "", "Please type a valid gene 2"),
            need(input$gene3 != "", "Please type a valid gene 3"),
            need(input$gene4 != "", "Please type a valid gene 4"),
            need(try(read.table(paste0(data_path, input$gene1,'_', input$gene2, file_suffix), sep='\t', header=FALSE)), 
                 paste0("There is no file for the pair ", input$gene1, ' and ', input$gene2)),
            need(try(read.table(paste0(data_path, input$gene3,'_', input$gene4, file_suffix), sep='\t', header=FALSE)), 
                 paste0("There is no file for the pair ", input$gene1, ' and ', input$gene2))
        )
        
        c(input$gene1, input$gene2, input$gene3, input$gene4)
    })
    
    output$boxplots <- renderPlot({
        genes = get_genes()
        file1 = paste0(data_path, genes[1],'_', genes[2], file_suffix)
        file2 = paste0(data_path, genes[3],'_', genes[4], file_suffix)
        data1 <- read.table(file1, sep='\t', header=FALSE)
        abs_data1 = log2(abs(data1[3]))
        
        data2 <- read.table(file2, sep='\t', header=FALSE)
        
        data1$pair = paste0(data1$V1,"_",data1$V2)
        data2$pair = paste0(data2$V1,"_",data2$V2)
        
        df_both = rbind(data1, data2)
        colnames(df_both)[3]="distance"
        df_both = df_both[,c(3,4)]
        df_both$distance=abs(df_both$distance)
        df_both$distance[mapply(is.infinite, df_both$distance)] = 0
        qplot(pair, distance, data = df_both, 
              geom=c("boxplot"), fill = pair,log="y")
    })
    
    
    output$densityplots <- renderPlot({
        genes = get_genes_density()
        file1 = paste0(data_path, genes[1],'_', genes[2], file_suffix)
        file2 = paste0(data_path, genes[3],'_', genes[4], file_suffix)
        data1 <- read.table(file1, sep='\t', header=FALSE)
        # abs_data1 = log2(abs(data1[3]))
        abs_data1 = abs(data1[3])
        
        data2 <- read.table(file2, sep='\t', header=FALSE)
        # abs_data2 = log2(abs(data2[3]))
        data1$pair = paste0(data1$V1,"_",data1$V2)
        data2$pair = paste0(data2$V1,"_",data2$V2)
        
        df_both = rbind(data1, data2)
        colnames(df_both)[3]="distance"
        df_both = df_both[,c(3,4)]
        df_both$distance=abs(df_both$distance)
        df_both$distance[mapply(is.infinite, df_both$distance)] = 0
        
        options(scipen=10000)
        
        qplot(distance, data = df_both, geom = "density", color = pair,y = ..scaled..,log="x")
    })
    
    output$stats1 <- renderPrint({
        genes = get_genes_ks()
        file1 = paste0(data_path, genes[1],'_', genes[2], file_suffix)
        data1 <- read.table(file1, sep='\t', header=FALSE)
        abs_data1 = log2(abs(data1[3]))
        
        # Replace -Inf caused by log-zero values
        abs_data1[mapply(is.infinite, abs_data1)] = 0
        
        paste0("log2 median intersummit distance between ", genes[1],' and ', genes[2], " is: ",median(unlist(abs_data1)))
        
    })
    output$stats2 <- renderPrint({
        genes = get_genes_ks()
        file2 = paste0(data_path, genes[3],'_', genes[4], file_suffix)
        
        data2 <- read.table(file2, sep='\t', header=FALSE)
        abs_data2 = log2(abs(data2[3]))
        
        # Replace -Inf caused by log-zero values
        abs_data2[mapply(is.infinite, abs_data2)] = 0
        
        paste0("log2 median intersummit distance between ", genes[3],' and ', genes[4], " is: ", median(unlist(abs_data2)))
        
    })
    
    
    
    
    output$ksTest <- renderPrint({
        genes = get_genes_ks()
        file1 = paste0(data_path, genes[1],'_', genes[2], file_suffix)
        file2 = paste0(data_path, genes[3],'_', genes[4], file_suffix)
        data1 <- read.table(file1, sep='\t', header=FALSE)
        abs_data1 = log2(abs(data1[3]))
        
        data2 <- read.table(file2, sep='\t', header=FALSE)
        abs_data2 = log2(abs(data2[3]))
        
        # Replace -Inf caused by log-zero values
        abs_data1[mapply(is.infinite, abs_data1)] = 0
        abs_data2[mapply(is.infinite, abs_data2)] = 0
        
        # ks = ks.test(unlist(data1[3]), unlist(data2[3]))
        ks = ks.test(unlist(abs_data1), unlist(abs_data2))
        
        paste0("median for data1 is: ",median(unlist(abs_data1)))
        paste0("median for data2 is: ",median(unlist(abs_data2)))
        
        
        
        ks
    })
}




# Run the application 
shinyApp(ui = ui, server = server)
