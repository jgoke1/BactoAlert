// Bacterial surveillance beacon; ASCII also works in plain workstation logs.
def startMessage(String pipelineVersion) {
    log.info """
    |  .-----------.
    |  |  (( ! ))  |   BACTOALERT
    |  '-----------'   Bacterial genomic surveillance
    |                  v${pipelineVersion}
    |""".stripMargin()
}
