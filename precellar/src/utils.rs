use std::ops::Range;

use std::sync::mpsc::{sync_channel, Receiver};
use std::thread::JoinHandle;
use noodles::fastq;
use regex::Regex;
use anyhow::{Result, anyhow};
use bstr::ByteSlice;

/// PrefetchIterator allows for prefetching items from an iterator into a buffer.
pub struct PrefetchIterator<T> {
    rx: Receiver<T>,
    // Keep the thread handle so we can join on Drop (nice for deterministic cleanup).
    handle: Option<JoinHandle<()>>,
}

impl<T: Send + 'static> PrefetchIterator<T> {
    pub fn new<I>(iter: I, buffer_size: usize) -> Self
    where
        I: IntoIterator<Item = T> + Send + 'static,
    {
        let (sender, receiver) = sync_channel(buffer_size);
        let handle = std::thread::spawn(move || {
            for item in iter {
                if sender.send(item).is_err() {
                    // If the receiver is dropped, we stop sending more items.
                    break;
                }
            }
        });
        Self { rx: receiver, handle: Some(handle) }
    }
}

impl<T: Send + 'static> Iterator for PrefetchIterator<T> {
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        self.rx.recv().ok()
    }
}

impl<T> Drop for PrefetchIterator<T> {
    fn drop(&mut self) {
        // Dropping rx will cause the producer to exit; join the thread to avoid detach.
        if let Some(handle) = self.handle.take() {
            let _ = handle.join();
        }
    }
}

pub struct PrefetchIterScoped<T> {
    rx: Receiver<T>,
}

impl<T> Iterator for PrefetchIterScoped<T> {
    type Item = T;
    fn next(&mut self) -> Option<T> {
        self.rx.recv().ok()
    }
}

/// Run `f` with a prefetching iterator, without requiring `'static` on `I` or `T`.
pub fn with_prefetch<'a, I, T, R, F>(iter: I, buffer: usize, f: F) -> R
where
    I: IntoIterator<Item = T> + Send + 'a,
    T: Send + 'a,
    F: FnOnce(PrefetchIterScoped<T>) -> R,
{
    std::thread::scope(|s| {
        let (tx, rx) = sync_channel::<T>(buffer);
        let _producer = s.spawn(move || {
            for item in iter {
                if tx.send(item).is_err() {
                    break;
                }
            }
        });
        let it = PrefetchIterScoped { rx };
        f(it) // consumes inside the scope; no `'static` needed
    })
}


#[derive(Debug, Clone)]
pub enum GroupIndex {
    Position(usize),
    Named(String),
}

impl From<usize> for GroupIndex {
    fn from(i: usize) -> Self {
        GroupIndex::Position(i)
    }
}

impl From<&str> for GroupIndex {
    fn from(s: &str) -> Self {
        GroupIndex::Named(s.to_string())
    }
}

impl From<String> for GroupIndex {
    fn from(s: String) -> Self {
        GroupIndex::Named(s)
    }
}

pub fn strip_fastq(
    mut fq: fastq::Record,
    regex: &Regex,
    group_indices: Option<&[GroupIndex]>,
    from_description: bool,
) -> Result<(fastq::Record, Option<Vec<String>>)> {
    let txt = if from_description {
        fq.description().to_str()?
    } else {
        fq.name().to_str()?
    };

    let (new_name, matches) = strip_from(txt, regex, group_indices)?;

    if from_description {
        *fq.description_mut() = new_name.into();
    } else {
        *fq.name_mut() = new_name.into();
    }

    Ok((fq, matches))
}

pub fn rev_compl_fastq_record(mut record: fastq::Record) -> fastq::Record {
    *record.quality_scores_mut() = record.quality_scores().iter().rev().copied().collect();
    *record.sequence_mut() = rev_compl(record.sequence());
    record
}

pub fn rev_compl(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&x| match x {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            _ => x,
        })
        .collect()
}

fn strip_from(
    txt: &str,
    regex: &Regex,
    group_indices: Option<&[GroupIndex]>,
) -> Result<(String, Option<Vec<String>>)> {
    let caps = regex.captures(txt).ok_or(anyhow!("No match found for: {}", txt))?;
    let new_name = remove_substrings(txt, caps.iter().skip(1).map(|m| m.unwrap().range()));
    let matches = if let Some(indices) = group_indices {
        let matches = indices.iter().map(|idx| match idx {
            GroupIndex::Position(i) => caps.get(*i+1).map(|m| m.as_str().to_string()),
            GroupIndex::Named(name) => caps.name(name).map(|m| m.as_str().to_string()),
        }).collect::<Option<Vec<_>>>().ok_or(anyhow!("Failed to extract barcodes"))?;
        Some(matches)
    } else {
        None
    };
    Ok((new_name, matches))
}

fn remove_substrings<I>(s: &str, ranges: I) -> String
where
    I: IntoIterator<Item = Range<usize>>,
{
    let mut result = String::with_capacity(s.len());
    let mut last_index = 0;

    for Range {start, end} in ranges.into_iter() {
        // Ensure valid ranges and skip invalid ones
        if start < end && end <= s.len() {
            result.push_str(&s[last_index..start]); // Append part before the range
            last_index = end; // Move past the removed range
        }
    }

    // Append remaining part after the last removed range
    result.push_str(&s[last_index..]);
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_strip_barcode() {
        let name = "A01535:24:HW2MMDSX2:2:1359:8513:3458:bd:69:Y6:10:TGATAGGTTG";
        let re = Regex::new(r"3458(:)(?<barcode>..:..:..:..)(:)(?<umi>[ATCG]+)$").unwrap();
        let indices = ["barcode".into(), "umi".into()];
        let (new_name, matches) = strip_from(name, &re, Some(indices.as_slice())).unwrap();
        let matches = matches.unwrap();
        assert_eq!(matches[0], "bd:69:Y6:10");
        assert_eq!(matches[1], "TGATAGGTTG");
        assert_eq!(new_name, "A01535:24:HW2MMDSX2:2:1359:8513:3458");
    }
}