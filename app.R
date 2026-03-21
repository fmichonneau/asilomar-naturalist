library(shiny)
library(bslib)
library(dplyr)
library(nanoparquet)
library(cachem)
library(httr2)
library(ggplot2)

# ââ Data prep (runs once at startup) âââââââââââââââââââââââââââââââââââââââââ

combined_data <- nanoparquet::read_parquet("combined_mry_data.parquet") |>
  mutate(observed_on = as.Date(observed_on))

taxa_table <- combined_data |>
  select(scientific_name, taxon_id, common_name, iconic_taxon_name) |>
  distinct()

# Top 100 per taxonomic group
top100 <- combined_data |>
  group_by(iconic_taxon_name, scientific_name) |>
  summarise(count = n(), .groups = "drop") |>
  slice_max(count, n = 100, by = iconic_taxon_name) |>
  left_join(taxa_table, by = c("scientific_name", "iconic_taxon_name")) |>
  select(scientific_name, common_name, iconic_taxon_name, taxon_id, count) |>
  arrange(iconic_taxon_name, desc(count))

# Friendly group label lookup
group_labels <- c(
  "Animalia" = "Invertebrates (Animalia)",
  "Mollusca" = "Snails, Slugs, Octopuses (Mollusca)",
  "Aves" = "Birds (Aves)",
  "Plantae" = "Plants (Plantae)",
  "Fungi" = "Fungi (Fungi)",
  "Actinopterygii" = "Ray-finned Fishes (Actinopterygii)",
  "Reptilia" = "Reptiles (Reptilia)",
  "Mammalia" = "Mammals (Mammalia)"
)

# All raw group values present in our top 100
raw_groups <- sort(unique(top100$iconic_taxon_name))

# Named vector: friendly label â raw value (for checkboxGroupInput)
choices_named <- setNames(
  raw_groups,
  vapply(
    raw_groups,
    function(g) {
      if (!is.na(group_labels[g])) {
        group_labels[[g]]
      } else {
        paste0(g, " (", g, ")")
      }
    },
    character(1)
  )
)

# ââ Image cache (disk, survives restarts) âââââââââââââââââââââââââââââââââââââ

img_cache <- cachem::cache_disk("image_cache")
wiki_cache <- cachem::cache_disk("wiki_cache")

fetch_taxon_info <- function(taxon_id) {
  cli::cli_alert_info(paste0("Fetching image for taxon_id ", taxon_id, " ..."))
  key <- as.character(taxon_id)

  img_cached <- img_cache$get(key)
  wiki_cached <- wiki_cache$get(key)

  if (!is.key_missing(img_cached) && !is.key_missing(wiki_cached)) {
    cli::cli_alert_success(paste0("Cache hit for taxon_id ", taxon_id))
    return(list(photo_url = img_cached, wikipedia_url = wiki_cached))
  }

  result <- tryCatch(
    {
      resp <- request(paste0(
        "https://api.inaturalist.org/v1/taxa/",
        taxon_id
      )) |>
        req_timeout(10) |>
        req_perform()
      body <- resp_body_json(resp)
      cli::cli_alert_success(paste0("Fetched image for taxon_id ", taxon_id))
      if (length(body$results) == 0) {
        return(list(photo_url = NA_character_, wikipedia_url = NA_character_))
      }
      taxon <- body$results[[1]]
      photo <- taxon$default_photo
      photo_url <- if (!is.null(photo$medium_url)) {
        photo$medium_url
      } else if (!is.null(photo$square_url)) {
        photo$square_url
      } else {
        NA_character_
      }
      # Fallback: if no taxon thumbnail, try the observations API
      if (is.na(photo_url)) {
        cli::cli_alert_info(paste0(
          "No taxon photo for ",
          taxon_id,
          ", trying observations..."
        ))
        obs_resp <- tryCatch(
          request("https://api.inaturalist.org/v1/observations") |>
            req_url_query(
              taxon_id = taxon_id,
              photos = "true",
              quality_grade = "research",
              order_by = "votes",
              per_page = 1L
            ) |>
            req_timeout(10) |>
            req_perform(),
          error = function(e) NULL
        )
        if (!is.null(obs_resp)) {
          obs_body <- resp_body_json(obs_resp)
          obs_photo <- tryCatch(
            obs_body$results[[1]]$photos[[1]]$url,
            error = function(e) NULL
          )
          if (!is.null(obs_photo)) {
            # iNat observation photo URLs use "square" size; upgrade to "medium"
            photo_url <- sub("/square\\.", "/medium.", obs_photo)
            cli::cli_alert_success(paste0(
              "Got observation photo for taxon_id ",
              taxon_id
            ))
          }
        }
      }

      wikipedia_url <- taxon$wikipedia_url %||% NA_character_
      list(photo_url = photo_url, wikipedia_url = wikipedia_url)
    },
    error = function(e) {
      list(photo_url = NA_character_, wikipedia_url = NA_character_)
    }
  )

  img_cache$set(key, result$photo_url)
  wiki_cache$set(key, result$wikipedia_url)
  result
}

fetch_taxon_photo <- function(taxon_id) {
  fetch_taxon_info(taxon_id)$photo_url
}

fetch_taxon_wiki <- function(taxon_id) {
  fetch_taxon_info(taxon_id)$wikipedia_url
}

# Pre-fetch all image URLs at startup
message("Fetching taxon images (cached after first run)...")
top100$photo_url <- vapply(top100$taxon_id, fetch_taxon_photo, character(1))
message("Images ready.")

# ââ UI ââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââ

ui <- page_sidebar(
  title = "Most Commonly Observed Species at Asilomar State Beach",
  theme = bs_theme(version = 5, bootswatch = "flatly"),

  sidebar = sidebar(
    width = 280,
    h6("Taxonomic Group"),
    tags$div(
      style = "margin-bottom: 6px; font-size: 0.85em;",
      actionLink("select_all", "Select all"),
      " / ",
      actionLink("select_none", "None")
    ),
    checkboxGroupInput(
      "groups",
      label = NULL,
      choices = choices_named,
      selected = raw_groups
    ),
    hr(),
    uiOutput("species_count")
  ),

  uiOutput("cards")
)

# ââ Server ââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââ

server <- function(input, output, session) {
  # Select all / none helpers
  observeEvent(input$select_all, {
    updateCheckboxGroupInput(session, "groups", selected = raw_groups)
  })
  observeEvent(input$select_none, {
    updateCheckboxGroupInput(session, "groups", selected = character(0))
  })

  # Reactive filtered data
  filtered <- reactive({
    req(input$groups)
    top100 |> filter(iconic_taxon_name %in% input$groups)
  })

  # Species count sub-header
  output$species_count <- renderUI({
    n <- nrow(filtered())
    tags$p(
      style = "color: #6c757d; font-size: 0.85em; margin: 0;",
      paste0(n, " species shown")
    )
  })

  # Ã¢Ã¢ Detail modal Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢

  # Reactive value holding the taxon row to display
  selected_species <- reactiveVal(NULL)

  observeEvent(
    input$show_detail,
    {
      req(input$show_detail)
      row <- top100[top100$taxon_id == input$show_detail, ][1, ]
      selected_species(row)

      obs <- combined_data |>
        filter(scientific_name == row$scientific_name) |>
        arrange(desc(observed_on))

      # Observation photo gallery (up to 6 non-NA image_urls)
      obs_photos <- obs |>
        filter(!is.na(image_url) & nchar(image_url) > 0) |>
        slice_head(n = 6) |>
        pull(image_url)

      gallery <- if (length(obs_photos) > 0) {
        photo_tags <- lapply(obs_photos, function(src) {
          tags$img(
            src = src,
            style = paste0(
              "width:150px; height:120px; object-fit:cover;",
              " border-radius:4px; flex-shrink:0;"
            ),
            alt = ""
          )
        })
        tags$div(
          style = paste0(
            "display:flex; flex-wrap:wrap; gap:6px;",
            " margin-bottom:16px;"
          ),
          !!!photo_tags
        )
      } else {
        tags$p(
          style = "color:#6c757d;",
          "No observation photos available."
        )
      }

      # Seasonal density plot
      today_yday <- as.integer(format(Sys.Date(), "%j"))
      obs_yday <- as.integer(format(obs$observed_on, "%j"))
      obs_yday <- obs_yday[!is.na(obs_yday)]

      # Effort-corrected relative frequency (week resolution)
      today_week <- as.integer(format(Sys.Date(), "%V"))

      effort <- combined_data |>
        mutate(week = as.integer(format(observed_on, "%V"))) |>
        count(week, name = "total_obs")

      rel_freq_df <- if (length(obs_yday) >= 2) {
        data.frame(
          week = as.integer(format(
            obs$observed_on[!is.na(obs$observed_on)],
            "%V"
          ))
        ) |>
          count(week, name = "n") |>
          left_join(effort, by = "week") |>
          mutate(relative_freq = n / total_obs)
      } else {
        NULL
      }

      density_plot <- if (length(obs_yday) >= 2) {
        renderPlot(
          {
            ggplot(data.frame(yday = obs_yday), aes(x = yday)) +
              geom_density(
                aes(y = after_stat(ndensity)),
                fill = "#2c7fb8",
                alpha = 0.4,
                color = "#2c7fb8"
              ) +
              geom_rug(alpha = 0.3, color = "#2c7fb8") +
              geom_vline(
                xintercept = today_yday,
                color = "#e63946",
                linewidth = 1,
                linetype = "dashed"
              ) +
              annotate(
                "text",
                x = today_yday,
                y = Inf,
                label = "Today",
                vjust = 1.5,
                hjust = if (today_yday > 300) 1.1 else -0.1,
                color = "#e63946",
                size = 3.5
              ) +
              scale_x_continuous(
                limits = c(1, 365),
                breaks = c(
                  1,
                  32,
                  60,
                  91,
                  121,
                  152,
                  182,
                  213,
                  244,
                  274,
                  305,
                  335
                ),
                labels = month.abb
              ) +
              labs(
                x = NULL,
                y = "Relative density",
                title = "Raw seasonal pattern"
              ) +
              theme_minimal(base_size = 12) +
              theme(panel.grid.minor = element_blank())
          },
          res = 120
        )
      } else {
        renderPlot({
          ggplot() +
            annotate(
              "text",
              x = 0.5,
              y = 0.5,
              label = "Not enough data for density plot"
            ) +
            theme_void()
        })
      }

      rel_freq_plot <- if (!is.null(rel_freq_df)) {
        renderPlot(
          {
            ggplot(rel_freq_df, aes(x = week, y = relative_freq)) +
              geom_col(fill = "#2c7fb8", alpha = 0.5, width = 0.8) +
              geom_smooth(
                se = FALSE,
                color = "#2c7fb8",
                linewidth = 1,
                method = "loess",
                formula = y ~ x
              ) +
              geom_vline(
                xintercept = today_week,
                color = "#e63946",
                linewidth = 1,
                linetype = "dashed"
              ) +
              annotate(
                "text",
                x = today_week,
                y = Inf,
                label = "Today",
                vjust = 1.5,
                hjust = if (today_week > 45) 1.1 else -0.1,
                color = "#e63946",
                size = 3.5
              ) +
              scale_x_continuous(
                limits = c(1, 53),
                breaks = c(1, 5, 9, 14, 18, 22, 27, 31, 35, 40, 44, 48),
                labels = month.abb
              ) +
              scale_y_continuous(
                labels = scales::percent_format(accuracy = 0.1)
              ) +
              labs(
                x = NULL,
                y = "% of all obs. that week",
                title = "Effort-corrected frequency"
              ) +
              theme_minimal(base_size = 12) +
              theme(panel.grid.minor = element_blank())
          },
          res = 120
        )
      } else {
        renderPlot({
          ggplot() +
            annotate(
              "text",
              x = 0.5,
              y = 0.5,
              label = "Not enough data"
            ) +
            theme_void()
        })
      }

      output$modal_density <- density_plot
      output$modal_rel_freq <- rel_freq_plot

      # Wikipedia link
      wiki_url <- fetch_taxon_wiki(row$taxon_id)

      wiki_link <- if (!is.na(wiki_url) && nchar(wiki_url) > 0) {
        tags$p(
          tags$a(
            href = wiki_url,
            target = "_blank",
            style = "font-size:0.9em;",
            paste0("Wikipedia: ", row$common_name %||% row$scientific_name)
          )
        )
      } else {
        NULL
      }

      # Observation list table
      obs_table <- obs |>
        select(observed_on, url) |>
        slice_head(n = 100)

      obs_rows <- lapply(seq_len(nrow(obs_table)), function(i) {
        tags$tr(
          tags$td(as.character(obs_table$observed_on[i])),
          tags$td(
            tags$a(
              href = obs_table$url[i],
              target = "_blank",
              style = "font-size:0.85em;",
              obs_table$url[i]
            )
          )
        )
      })

      inat_url <- paste0("https://www.inaturalist.org/taxa/", row$taxon_id)
      modal_title <- tags$span(
        tags$a(
          href = inat_url,
          target = "_blank",
          style = "color:inherit; text-decoration:none;",
          row$common_name %||% row$scientific_name
        ),
        tags$small(
          style = "color:#6c757d; font-style:italic; margin-left:8px;",
          row$scientific_name
        )
      )

      showModal(modalDialog(
        title = modal_title,
        size = "xl",
        easyClose = TRUE,
        footer = modalButton("Close"),
        gallery,
        wiki_link,
        h6("Seasonal pattern"),
        tags$div(
          style = "display:flex; gap:12px;",
          tags$div(
            style = "flex:1;",
            plotOutput("modal_density", height = "200px")
          ),
          tags$div(
            style = "flex:1;",
            plotOutput("modal_rel_freq", height = "200px")
          )
        ),
        h6(
          style = "margin-top:16px;",
          paste0("Observations (", format(row$count, big.mark = ","), " total)")
        ),
        tags$div(
          style = "max-height:300px; overflow-y:auto;",
          tags$table(
            class = "table table-sm table-striped",
            tags$thead(
              tags$tr(tags$th("Date"), tags$th("iNaturalist URL"))
            ),
            tags$tbody(!!!obs_rows)
          )
        )
      ))
    },
    ignoreNULL = TRUE,
    ignoreInit = TRUE
  )

  # Ã¢Ã¢ Card grid Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢Ã¢

  output$cards <- renderUI({
    df <- filtered()
    if (nrow(df) == 0) {
      return(p("No species match the selected filters."))
    }

    cards <- lapply(seq_len(nrow(df)), function(i) {
      row <- df[i, ]
      inat_url <- paste0("https://www.inaturalist.org/taxa/", row$taxon_id)
      photo_url <- row$photo_url
      count_fmt <- format(row$count, big.mark = ",")

      # Thumbnail or placeholder
      img_tag <- if (!is.na(photo_url)) {
        tags$a(
          href = inat_url,
          target = "_blank",
          tags$img(
            src = photo_url,
            style = "width:100%; height:150px; object-fit:cover; border-radius:4px 4px 0 0;",
            alt = row$common_name
          )
        )
      } else {
        tags$div(
          style = paste0(
            "width:100%; height:150px; background:#e9ecef; border-radius:4px 4px 0 0;",
            "display:flex; align-items:center; justify-content:center; color:#adb5bd;"
          ),
          "No image"
        )
      }

      card(
        full_screen = FALSE,
        style = "padding: 0; overflow: hidden;",
        img_tag,
        card_body(
          style = "padding: 10px;",
          tags$a(
            href = inat_url,
            target = "_blank",
            style = "font-weight:600; color: inherit; text-decoration: none;",
            row$common_name %||% row$scientific_name
          ),
          tags$br(),
          tags$em(
            style = "color:#6c757d; font-size:0.85em;",
            row$scientific_name
          ),
          tags$br(),
          tags$span(
            class = "badge bg-secondary mt-1",
            style = "cursor:pointer;",
            onclick = sprintf(
              "Shiny.setInputValue('show_detail', %s, {priority: 'event'})",
              row$taxon_id
            ),
            paste0(count_fmt, " observations")
          )
        )
      )
    })

    layout_column_wrap(width = "200px", !!!cards)
  })
}

shinyApp(ui, server)
