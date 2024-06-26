% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ExtractImageFeature.R
\name{ExtractFeatures}
\alias{ExtractFeatures}
\title{extract image features}
\usage{
ExtractFeatures(
  imgFile,
  positionTable,
  patchSize = NULL,
  scaleFactor = 1,
  saveRawImage = FALSE,
  savePatchOnImage = NULL
)
}
\arguments{
\item{imgFile}{the histology image file}

\item{positionTable}{spots coordinate information}

\item{patchSize}{the diameter of patches. Suggest 256 for ST platform.
For Visium platform, the spot diameter in \code{scalefactors_json.json} is
suggested.}

\item{scaleFactor}{A scaling factor that converts pixel positions in the
original, full-resolution image to pixel positions in the image file used
here. 1 if used the raw image.}

\item{saveRawImage}{logical value to define wether save raw image to file.}

\item{savePatchOnImage}{path and file name of plotting patches on the image.}
}
\value{
A list of two matrices
}
\description{
extract image features
}
\examples{
\dontrun{
# ST platform
imageFeatures.list <-
  ExtractFeatures(
    imgFile = rawimg,
    positionTable = positionTable,
    patchSize = patchsize,
    rawImage = T,
    savePatchOnImage = paste0("sampleName_rawpatch_", patchsize, ".png")
  )
imageFeatures <- imageFeatures.list[["ImageFeature"]]
RGBquantile <- imageFeatures.list[["RGBquantile"]]
}

}
