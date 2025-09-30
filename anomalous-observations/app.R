library(shiny)
library(bslib)

ui <- page_navbar(
  sidebar = sidebar(
    selectInput(inputId = "phenophase",
                label = "Phenophase class",
                choices = c(1, 3, 6, 7))
  ),
  title = "Anomalous observations",
  bg = "#2D89C8",
  inverse = TRUE,
  nav_spacer(), 
  nav_panel(title = "Map", p("Leaflet map with locations")),
  nav_panel(title = "Table", p("Details about anomalous observations"))
)

server <- function(input, output) {}

# Run the application
shinyApp(ui = ui, server = server)