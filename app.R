#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
# library (shinydashboard)

library(dplyr)
library(ggplot2)
# Define UI for application
ui <- fluidPage(

  # Application title
  titlePanel("Transcription factor visualization in gene promoters"),

  textOutput('welcome_text'),
  # load files
  sidebarLayout(
    sidebarPanel(
      textOutput('text_for_files'),
  fileInput("upload", "Upload a file", accept = c("text/csv",
                                                  "text/comma-separated-values,text/plain",
                                                  ".csv")),
  
  downloadButton("download", label = "Download sample data"),
  actionButton('sample_run', 'Run using sample data'),
  textInput('gene_name', 'Enter gene name for plot title (optional)'),
  numericInput('prom_length', 'Enter promoter length (default 1500)', 
               value = 1500, max = 2500, min = 0),
  tableOutput("head"),
  fluidRow(
    
    actionButton("result1","Generate Result"),
    downloadButton('download_plot','Download Plot')
    ,plotOutput("plot_sample", height = '550px', width = '710px')
  ),
  'Panel for sample plot'
  ),
  mainPanel(
    textOutput('main_header'),
    # plotOutput("tf_plot", height = paste0(1500 - 1500/2, 'px'), width = "1200px")
    uiOutput("ui_plot")
  )
)
)

# height = 1500-prom_length/2, width = 1500-prom_length, units = 'px'

# Define server logic required to draw a TF ggplot2
server <- function(input, output) {
  output$welcome_text <- renderText(('It is alpha-version'))
  output$text_for_files <- renderText('Files could be tab- or comma-separated in .txt/.csv format without headers')
  output$main_header <- renderText('Your plot')
  output$download <- downloadHandler(
    filename <- function() {
      "template_output.txt"
    },
    content <- function(file) {
      file.copy("./data/template.txt", file)
    }
  )
  df <- reactive({
    req(input$upload)
    
    ext <- tools::file_ext(input$upload$name)
    switch(ext,
           csv = vroom::vroom(input$upload$datapath, delim = ",", col_names = c('TF', 'rel_beg', 'rel_end')),
           tsv = vroom::vroom(input$upload$datapath, delim = "\t", col_names = c('TF', 'rel_beg', 'rel_end')),
           txt = vroom::vroom(input$upload$datapath, delim = "\t", col_names = c('TF', 'rel_beg', 'rel_end')),
           validate("Invalid file; Please upload a .csv/.tsv/.txt file")
    )
  })
  
  output$head <- renderTable({
    head(df(), 10)
  })

  x_breaks_prom <- reactive(c(seq(-input$prom_length-500,400,300)[seq(-input$prom_length-500,400,300)!= 0]))
  
  # 
  data <- reactiveValues()
  observeEvent(input$sample_run, {
    sample_df <- reactive(vroom::vroom('data/template.txt', delim = "\t", col_names = c('TF', 'rel_beg', 'rel_end')))
    data$sample_plot <- sample_df() %>%
      # mutate(height = 1+(1:nrow(df)))%>%
      group_by(TF) %>% # to group_by
      mutate(height = cur_group_id()) %>% # to create column with same height for same TF
      ggplot(aes(x = height, group = TF, color = TF, fill = TF))+
      geom_rect(mapping=aes(xmin=rel_beg, xmax=rel_end,
                            ymin=height, ymax=height-0.2), alpha=0.5, show.legend = FALSE)+ # draw TF directly
      
      geom_text(aes(rel_beg-30, height-0.08, label = rel_beg), size = 2.5, color = 'black', show.legend = FALSE)+ # draw coordinate of beginning
      geom_text(aes(rel_end+30, height-0.08, label = rel_end), size = 2.5, color = 'black', show.legend = FALSE)+ # draw coordinate of ending
      geom_text(aes((rel_end+rel_beg)/2, height+0.5, label = TF), color = 'black',
                size = 2.5, show.legend = FALSE)+ #draw TF label
      scale_x_continuous(name = 'promoter', limits = c(-1500-500,350),
                         breaks = c(x_breaks_prom(),1),
                         labels = c(x_breaks_prom(),'+1'),
                         # sec.axis = dup_axis( name = paste0(c('genome coordinates',chr,':'),collapse = ' '),
                         #                      labels = x_breaks_gc)
      )+
      
      scale_y_continuous(limits = c(-0.5,(nrow(sample_df())+1+1)), expand = c(0,0))+
      theme_classic()+
      ggtitle('Sample gene')+
      # draw arrow as TSS
      geom_segment(aes(x = 1, y = -0.5, xend = 1, yend = 1), color = 'black')+
      geom_segment(aes(x = 1, y = 1, xend = 100, yend = 1), arrow = arrow(length = unit(0.02,units = "npc")),
                   size = 0.8, color = 'black')+
      geom_text(aes(20, 2.1, label = 'TSS'), color = 'black', size = 3)+
      theme(plot.title = element_text(hjust = 0.5, face = 'italic'),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.line.y = element_blank())
    # plot_sample  <- renderPlot({  data$sample_plot })
  })
  output$plot_sample <- renderPlot({data$sample_plot}, res = 96)
  
  observeEvent(input$result1,{
    data$plot <- df() %>%
      # mutate(height = 1+(1:nrow(df)))%>%
      group_by(TF) %>% # to group_by
      mutate(height = cur_group_id()) %>% # to create column with same height for same TF
      ggplot(aes(x = height, group = TF, color = TF, fill = TF))+
      geom_rect(mapping=aes(xmin=rel_beg, xmax=rel_end,
                            ymin=height, ymax=height-0.2), alpha=0.5, show.legend = FALSE)+ # draw TF directly

      geom_text(aes(rel_beg-30, height-0.08, label = rel_beg), size = 2.5, color = 'black', show.legend = FALSE)+ # draw coordinate of beginning
      geom_text(aes(rel_end+30, height-0.08, label = rel_end), size = 2.5, color = 'black', show.legend = FALSE)+ # draw coordinate of ending
      geom_text(aes((rel_end+rel_beg)/2, height+0.5, label = TF), color = 'black',
                size = 2.5, show.legend = FALSE)+ #draw TF label
      scale_x_continuous(name = 'promoter', limits = c(-input$prom_length-500,350),
                         breaks = c(x_breaks_prom(),1),
                         labels = c(x_breaks_prom(),'+1'),
                         # sec.axis = dup_axis( name = paste0(c('genome coordinates',chr,':'),collapse = ' '),
                         #                      labels = x_breaks_gc)
      )+

      scale_y_continuous(limits = c(-0.5,(nrow(df())+1+1)), expand = c(0,0))+
      theme_classic()+
      ggtitle(input$gene_name)+
      # draw arrow as TSS
      geom_segment(aes(x = 1, y = -0.5, xend = 1, yend = 1), color = 'black')+
      geom_segment(aes(x = 1, y = 1, xend = 100, yend = 1), arrow = arrow(length = unit(0.02,units = "npc")),
                   size = 0.8, color = 'black')+
      geom_text(aes(20, 2.1, label = 'TSS'), color = 'black', size = 3)+
      theme(plot.title = element_text(hjust = 0.5, face = 'italic'),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.line.y = element_blank())
    })

  plot3  <- renderPlot({  data$plot })
  output$tf_plot <- renderPlot({data$plot}, res = 96)
  
  output$ui_plot <- renderUI({
    plotOutput("tf_plot", height = (1500+input$prom_length/2)/3, width = (1500+input$prom_length)/3)
  })
  
  output$download_plot <- downloadHandler(
    filename = function(){paste("TF_plot",'.png',sep='')},
    content = function(file){
      ggsave(file,plot=data$plot, height = (1500+input$prom_length/2), width = (1500+input$prom_length), dpi = 300, units = 'px')
    }
  )
  # 


}

# Run the application 
shinyApp(ui = ui, server = server)
