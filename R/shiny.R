#' @export
shiny_ui <- function()
{
  fluidPage(
    fileInput("upload", NULL, buttonLabel = "Upload...", multiple = TRUE),
    DT::DTOutput("files"),
    tableOutput("selected_file_table")
  )
}

#' @export
shiny_server <- function(input, output, session) {
  output$files <- DT::renderDT({
    DT::datatable(input$upload, selection = c("single"))
  })

  # read all the uploaded files
  all_files <- reactive({
    req(input$upload)
    purrr::map(input$upload$datapath, readr::read_csv) %>%
      purrr::set_names(input$upload$name)
  })

  #select a row in DT files and display the corresponding table
  output$selected_file_table <- renderTable({
    req(input$upload)
    req(input$files_rows_selected)

    all_files()[[
    input$upload$name[[input$files_rows_selected]]
    ]]
  })
}

#shinyApp(ui, server)
