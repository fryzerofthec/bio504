options(timeout = 300)

accession <- "NG_000006.1"
genes     <- c("HBA1", "HBA2", "HBM", "HBQ1", "HBZ")
outdir    <- "fasta_globines"
zipname   <- "alpha_globines_fasta.zip"
line_w    <- 70

# ---------- utilitaires ----------
efetch <- function(db, id, rettype) {
  u <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=", db,
              "&id=", id, "&rettype=", rettype, "&retmode=text")
  con <- url(u)
  on.exit(close(con))
  readLines(con, warn = FALSE)
}

wrap_seq <- function(s, w = line_w) {
  starts <- seq(1, nchar(s), by = w)
  substring(s, starts, starts + w - 1)
}

write_fasta <- function(path, header, s) {
  writeLines(c(paste0(">", header), wrap_seq(s)), path)
}

revcomp <- function(x) {
  x <- chartr("ACGTRYMKBDHVN", "TGCAYRKMVHDBN", x)
  paste(rev(strsplit(x, "")[[1]]), collapse = "")
}

translate <- function(cds) {
  b  <- c("T", "C", "A", "G")
  aa <- strsplit("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "")[[1]]
  n  <- nchar(cds) %/% 3
  st <- (0:(n - 1)) * 3 + 1
  codons <- substring(cds, st, st + 2)
  res <- vapply(codons, function(cd) {
    x <- match(strsplit(cd, "")[[1]], b) - 1
    if (anyNA(x)) "X" else aa[x[1] * 16 + x[2] * 4 + x[3] + 1]
  }, character(1))
  res <- unname(res)
  stop_pos <- match("*", res)
  if (!is.na(stop_pos)) res <- res[seq_len(stop_pos - 1)]
  paste(res, collapse = "")
}

# ---------- 1. séquence complète----------
message("Téléchargement de la séquence ", accession, " ...")
fa  <- efetch("nuccore", accession, "fasta")
acc <- strsplit(sub("^>", "", fa[1]), " ")[[1]][1]          
nuc <- toupper(gsub("\\s", "", paste(fa[-1], collapse = "")))
message("  ", acc, " : ", nchar(nuc), " pb")

# ---------- 2. table des annotations (feature table) ----------
message("Téléchargement des annotations ...")
ft <- efetch("nuccore", accession, "ft")

parse_ft <- function(lines) {
  feats <- list(); n <- 0
  for (ln in lines) {
    if (!nzchar(ln) || startsWith(ln, ">Feature")) next
    f <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    if (startsWith(ln, "\t\t\t")) {                         
      if (n == 0) next
      key <- f[4]
      val <- if (length(f) >= 5) f[5] else ""
      feats[[n]]$q[[key]] <- c(feats[[n]]$q[[key]], val)
    } else {                                                
      s <- as.numeric(gsub("[<>]", "", f[1]))
      e <- as.numeric(gsub("[<>]", "", f[2]))
      type <- if (length(f) >= 3 && nzchar(f[3])) f[3] else NA
      if (!is.na(type)) {
        n <- n + 1
        feats[[n]] <- list(type = type, iv = list(c(s, e)), q = list())
      } else if (n > 0) {
        feats[[n]]$iv[[length(feats[[n]]$iv) + 1]] <- c(s, e)
      }
    }
  }
  feats
}
feats <- parse_ft(ft)
message("  ", length(feats), " annotations lues")

extract_seq <- function(feat) {
  parts <- vapply(feat$iv, function(v) {
    s <- substring(nuc, min(v), max(v))
    if (v[1] > v[2]) revcomp(s) else s                      
  }, character(1))
  paste(parts, collapse = "")
}

loc_txt <- function(feat) {
  all_pos <- unlist(feat$iv)
  strand  <- if (feat$iv[[1]][1] > feat$iv[[1]][2]) "minus" else "plus"
  paste0(acc, ":", min(all_pos), "-", max(all_pos), " (", strand, " strand)")
}

pick <- function(type, g) {
  Filter(function(f) f$type == type && g %in% f$q[["gene"]], feats)
}

# ---------- 3. génération des fichiers ----------
dir.create(outdir, showWarnings = FALSE)
written <- character(0)

for (g in genes) {
  # gène
  gf <- pick("gene", g)
  if (length(gf) == 0) stop("Gène ", g, " introuvable dans ", acc)
  gf <- gf[[1]]
  p <- file.path(outdir, paste0(g, "_gene.txt"))
  write_fasta(p, paste0(loc_txt(gf), " ", g, " gene, genomic sequence"), extract_seq(gf))
  written <- c(written, p)

  # ARNm
  mf <- Filter(function(f) !is.null(f$q[["transcript_id"]]), pick("mRNA", g))
  if (length(mf) == 0) stop("ARNm de ", g, " introuvable")
  if (length(mf) > 1) message("  [!] ", g, " : ", length(mf), " ARNm, le premier est utilisé")
  mf  <- mf[[1]]
  tid <- mf$q[["transcript_id"]][1]
  p <- file.path(outdir, paste0(g, "_mRNA.txt"))
  write_fasta(p, paste0(tid, " ", g, " mRNA (DNA alphabet, T not replaced by U) [", loc_txt(mf), "]"),
              extract_seq(mf))
  written <- c(written, p)

  # CDS
  cf <- pick("CDS", g)
  if (length(cf) == 0) stop("CDS de ", g, " introuvable")
  cf  <- cf[[1]]
  cds <- extract_seq(cf)
  p <- file.path(outdir, paste0(g, "_CDS.txt"))
  write_fasta(p, paste0(loc_txt(cf), " ", g, " CDS"), cds)
  written <- c(written, p)

  # protéine (traduction de la CDS)
  pid  <- if (!is.null(cf$q[["protein_id"]])) cf$q[["protein_id"]][1] else "protein"
  prot <- translate(cds)
  p <- file.path(outdir, paste0(g, "_protein.txt"))
  write_fasta(p, paste0(pid, " ", g, " protein"), prot)
  written <- c(written, p)

  # contrôle : comparaison avec la protéine RefSeq de NCBI
  if (pid != "protein") {
    Sys.sleep(0.4)
    ref <- tryCatch(paste(efetch("protein", pid, "fasta")[-1], collapse = ""),
                    error = function(e) NA_character_)
    if (!is.na(ref) && ref != prot) {
      warning(g, " : la traduction de la CDS diffère de la protéine NCBI ", pid)
    }
  }
  message("  ", g, " : OK")
}

# ---------- 4. zip ----------
ok <- tryCatch({
  if (!requireNamespace("zip", quietly = TRUE)) install.packages("zip")
  zip::zipr(zipname, files = written)
  TRUE
}, error = function(e) FALSE)
if (!ok) ok <- tryCatch({ utils::zip(zipname, written); TRUE }, error = function(e) FALSE)

message(length(written), " fichiers écrits dans '", file.path(getwd(), outdir), "'")
if (ok) message("Archive : ", file.path(getwd(), zipname)) else
  message("Zip non créé : compressez manuellement le dossier '", outdir, "'.")
