use log::warn;
use noodles::fastq;
use regex::Regex;
use anyhow::{Result, anyhow};
use bstr::ByteSlice;

pub fn strip_barcode_from_read_name(
    mut fq: fastq::Record,
    regex: &Regex,
    left_add: usize,
    right_add: usize,
) -> Result<(fastq::Record, String)> {
    let read_name = fq.name().to_str()?;
    let (barcode, name) = remove_barcode(read_name, regex, left_add, right_add)?;
    *fq.name_mut() = name.into();
    Ok((fq, barcode))
}

fn remove_barcode(name: &str, re: &Regex, left_add: usize, right_add: usize) -> Result<(String, String)> {
    let mut mat = re.captures(name)
        .and_then(|x| x.get(1))
        .ok_or(anyhow!("The regex must contain exactly one capturing group matching the barcode"))?
        .range();
    let barcode = name.get(mat.clone()).unwrap().to_string();
    if barcode.is_empty() {
        warn!("regex match is empty for read name: {}", name);
        Ok((barcode, name.to_string()))
    } else {
        mat = (mat.start - left_add)..(mat.end + right_add);
        let new_name = name.get(..mat.start).unwrap().to_string() + name.get(mat.end..).unwrap();
        Ok((barcode, new_name))
    }

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_strip_barcode() {
        let name = "A01535:24:HW2MMDSX2:2:1359:8513:3458:bd:69:Y6:10:TGATAGGTTG";
        let re = Regex::new(r"(..:..:..:..):\w+$").unwrap();
        let (barcode, new_name) = remove_barcode(name, &re, 1, 1).unwrap();
        assert_eq!(barcode, "bd:69:Y6:10");
        assert_eq!(new_name, "A01535:24:HW2MMDSX2:2:1359:8513:3458TGATAGGTTG");
    }
}