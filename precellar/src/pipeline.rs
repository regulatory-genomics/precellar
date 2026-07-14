//! Reusable stages for annotated FASTQ batches.

use crate::align::AnnotatedFastq;
use anyhow::Result;
use serde_json::Value;

pub type AnnotatedFastqBatch = Vec<AnnotatedFastq>;

#[derive(Debug)]
pub struct MiddlewareQcReport {
    pub name: String,
    pub metrics: Value,
}

/// A stateful transformation applied to annotated FASTQ batches.
pub trait FastqStage: Send {
    fn process(&mut self, batch: AnnotatedFastqBatch) -> Result<AnnotatedFastqBatch>;

    fn finish(&mut self) -> Result<()> {
        Ok(())
    }

    fn report(&self) -> Option<MiddlewareQcReport> {
        None
    }
}

/// Ordered middleware stages applied to annotated FASTQ batches.
#[derive(Default)]
pub struct FastqStagePipeline {
    stages: Vec<Box<dyn FastqStage>>,
    finished: bool,
}

impl FastqStagePipeline {
    pub fn push_stage<S>(&mut self, stage: S)
    where
        S: FastqStage + 'static,
    {
        self.stages.push(Box::new(stage));
    }

    pub fn process(&mut self, mut batch: AnnotatedFastqBatch) -> Result<AnnotatedFastqBatch> {
        for stage in &mut self.stages {
            batch = stage.process(batch)?;
        }
        Ok(batch)
    }

    pub fn finish(&mut self) -> Result<()> {
        if self.finished {
            return Ok(());
        }
        for stage in &mut self.stages {
            stage.finish()?;
        }
        self.finished = true;
        Ok(())
    }

    pub fn reports(&self) -> Vec<MiddlewareQcReport> {
        self.stages
            .iter()
            .filter_map(|stage| stage.report())
            .collect()
    }
}
