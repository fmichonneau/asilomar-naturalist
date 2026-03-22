library(shiny)
library(bslib)
library(dplyr)
library(nanoparquet)
library(cachem)
library(httr2)
library(ggplot2)


# ââ Data prep (runs once at startup) âââââââââââââââââââââââââââââââââââââââââ

combined_data <- nanoparquet::read_parquet("combined_mry_data.parquet") |>
  mutate(observed_on = as.Date(observed_on)) |>
  filter(nzchar(iconic_taxon_name))

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

# Count taxa per group, ordered descending
group_counts <- top100 |>
  count(iconic_taxon_name, name = "n_taxa") |>
  arrange(desc(n_taxa))

# Named vector: friendly label -> raw value (for checkboxGroupInput)
# ordered by descending taxa count, with count in label
choices_named <- setNames(
  group_counts$iconic_taxon_name,
  vapply(
    group_counts$iconic_taxon_name,
    function(g) {
      n <- group_counts$n_taxa[group_counts$iconic_taxon_name == g]
      label <- if (!is.na(group_labels[g])) group_labels[[g]] else g
      paste0(label, " (", n, ")")
    },
    character(1)
  )
)
raw_groups <- unname(choices_named)

# ââ Image cache (disk, survives restarts) âââââââââââââââââââââââââââââââââââââ

img_dir <- "image_cache_files"
dir.create(img_dir, showWarnings = FALSE)
addResourcePath("taxon_images", img_dir)

# Stores iNaturalist photo URLs; expires after 90 days to refresh from API
img_url_cache <- cachem::cache_disk("img_cache_url", max_age = 90 * 24 * 3600)
wiki_cache <- cachem::cache_disk("wiki_cache", max_age = 90 * 24 * 3600)

fetch_taxon_info <- function(taxon_id) {
  cli::cli_alert_info(paste0("Fetching image for taxon_id ", taxon_id, " ..."))
  key <- as.character(taxon_id)
  img_path <- file.path(img_dir, paste0(key, ".jpg"))

  cached_url <- img_url_cache$get(key)
  cached_wiki <- wiki_cache$get(key)

  # Full cache hit: URL not expired AND file already downloaded
  if (
    !cachem::is.key_missing(cached_url) &&
      file.exists(img_path) &&
      !cachem::is.key_missing(cached_wiki)
  ) {
    cli::cli_alert_success(paste0("Cache hit for taxon_id ", taxon_id))
    return(list(
      photo_url = paste0("taxon_images/", key, ".jpg"),
      wikipedia_url = cached_wiki
    ))
  }

  # Fetch fresh data from iNaturalist API
  result <- tryCatch(
    {
      resp <- request(paste0(
        "https://api.inaturalist.org/v1/taxa/",
        taxon_id
      )) |>
        req_throttle(capacity = 50, fill_time_s = 60) |>
        req_retry(
          max_tries = 5,
          is_transient = \(resp) resp_status(resp) %in% c(429, 503),
          backoff = \(i) min(2^i, 60)
        ) |>
        req_perform()
      body <- resp_body_json(resp)
      cli::cli_alert_success(paste0(
        "Fetched taxon info for taxon_id ",
        taxon_id
      ))
      if (length(body$results) == 0) {
        return(list(
          photo_url = NA_character_,
          wikipedia_url = NA_character_
        ))
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
      list(
        photo_url = photo_url,
        wikipedia_url = taxon$wikipedia_url %||% NA_character_
      )
    },
    error = function(e) {
      cli::cli_alert_warning(paste0(
        "Error fetching taxon info for taxon_id ",
        taxon_id,
        ": ",
        conditionMessage(e)
      ))
      list(photo_url = NA_character_, wikipedia_url = NA_character_)
    }
  )

  # Download image to local cache
  local_url <- if (!is.na(result$photo_url)) {
    tryCatch(
      {
        download.file(result$photo_url, img_path, quiet = TRUE, mode = "wb")
        paste0("taxon_images/", key, ".jpg")
      },
      error = function(e) {
        cli::cli_alert_warning(paste0(
          "Failed to download image for taxon_id ",
          taxon_id,
          ": ",
          e$message
        ))
        NA_character_
      }
    )
  } else {
    cli::cli_alert_warning(paste0("No image found for taxon_id ", taxon_id))
    NA_character_
  }

  # Cache the URL and wiki link for next time
  img_url_cache$set(key, result$photo_url)
  wiki_cache$set(key, result$wikipedia_url)

  list(photo_url = local_url, wikipedia_url = result$wikipedia_url)
}

# Pre-fetch all image URLs at startup
message("Fetching taxon images (cached after first run)...")
top100$photo_url <- vapply(
  top100$taxon_id,
  function(id) fetch_taxon_info(id)$photo_url,
  character(1)
)
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
      obs_with_photos <- obs |>
        filter(!is.na(image_url) & nchar(image_url) > 0) |>
        slice_head(n = 12)

      gallery <- if (nrow(obs_with_photos) > 0) {
        photo_tags <- lapply(seq_len(nrow(obs_with_photos)), function(j) {
          src <- obs_with_photos$image_url[j]
          obs_url <- obs_with_photos$url[j]
          img <- tags$img(
            src = src,
            style = paste0(
              "width:150px; height:120px; object-fit:cover;",
              " border-radius:4px; flex-shrink:0;"
            ),
            alt = ""
          )
          if (!is.na(obs_url) && nchar(obs_url) > 0) {
            tags$a(
              href = obs_url,
              target = "_blank",
              title = "View observation on iNaturalist",
              style = "display:inline-block; flex-shrink:0;",
              img
            )
          } else {
            img
          }
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
      wiki_url <- fetch_taxon_info(row$taxon_id)$wikipedia_url

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

    # Preserve sidebar order: use choices_named values (which are raw group names)
    # filtered to only the selected groups, in that order
    group_order <- choices_named[choices_named %in% df$iconic_taxon_name]

    sections <- lapply(group_order, function(grp) {
      grp_df <- df[df$iconic_taxon_name == grp, ]
      label <- if (!is.na(group_labels[grp])) group_labels[[grp]] else grp

      cards <- lapply(seq_len(nrow(grp_df)), function(i) {
        row <- grp_df[i, ]
        inat_url <- paste0("https://www.inaturalist.org/taxa/", row$taxon_id)
        photo_url <- row$photo_url
        count_fmt <- format(row$count, big.mark = ",")

        show_detail_js <- sprintf(
          "Shiny.setInputValue('show_detail', %s, {priority: 'event'})",
          row$taxon_id
        )

        img_tag <- if (!is.na(photo_url)) {
          tags$div(
            style = "cursor:pointer;",
            onclick = show_detail_js,
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
              "display:flex; align-items:center; justify-content:center; color:#adb5bd;",
              "cursor:pointer;"
            ),
            onclick = show_detail_js,
            "No image"
          )
        }

        card(
          full_screen = FALSE,
          style = "padding: 0; overflow: hidden;",
          img_tag,
          card_body(
            style = "padding: 10px;",
            tags$span(
              style = "font-weight:600; cursor:pointer;",
              onclick = show_detail_js,
              row$common_name %||% row$scientific_name
            ),
            tags$br(),
            tags$em(
              style = "color:#6c757d; font-size:0.85em;",
              row$scientific_name
            ),
            tags$a(
              href = inat_url,
              target = "_blank",
              style = "font-size:0.78em; color:#6c757d; text-decoration:none;",
              "View on iNaturalist \u2197"
            ),
            tags$br(),
            tags$span(
              class = "badge bg-secondary mt-1",
              style = "cursor:pointer;",
              onclick = show_detail_js,
              paste0(count_fmt, " observations")
            )
          )
        )
      })

      tagList(
        h4(label, style = "margin-top: 1.2em; margin-bottom: 0.5em;"),
        layout_column_wrap(width = "200px", !!!cards)
      )
    })

    tagList(!!!sections)
  })
}

shinyApp(ui, server)
