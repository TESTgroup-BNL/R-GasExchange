\name{arrhenius.scaling}
\alias{arrhenius.scaling}
\title{An Arrhenius temperature scaling function for temperature dependent parameters}
\usage{
  arrhenius.scaling(temp1 = NULL, temp2 = NULL, E = NULL,
    Val.temp1 = NULL)
}
\arguments{
  \item{temp1}{the temperature (in degrees C) of the
  observation at the base temperature}

  \item{temp2}{the desired temperature (in degrees C) to
  scale the observation}

  \item{E}{activation energy of the Arrhenius equation}

  \item{Val.temp1}{the value of the observation at the base
  temperature}
}
\description{
  Arrhenius temperature scaling function
}
\examples{
arrhenius.scaling(temp1=15,temp2=25,E=54.08,Val.temp1=40)
}
\author{
  Shawn P. Serbin
}

