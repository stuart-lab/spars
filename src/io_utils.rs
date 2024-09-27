use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};

pub fn open_read(filename: &str) -> std::io::Result<Box<dyn Read>> {
    let file = File::open(filename)?;
    if filename.ends_with(".gz") {
        let decoder = GzDecoder::new(file);
        Ok(Box::new(decoder))
    } else {
        Ok(Box::new(file))
    }
}

pub fn open_write(filename: &str) -> std::io::Result<Box<dyn Write>> {
    let file = File::create(filename)?;
    if filename.ends_with(".gz") {
        let encoder = GzEncoder::new(file, Compression::default());
        Ok(Box::new(encoder))
    } else {
        Ok(Box::new(file))
    }
}

pub fn get_reader(filename: &str) -> std::io::Result<Box<dyn BufRead>> {
    let reader = BufReader::new(open_read(filename)?);
    Ok(Box::new(reader))
}

pub fn get_writer(filename: &str) -> std::io::Result<Box<dyn Write>> {
    let writer = BufWriter::new(open_write(filename)?);
    Ok(Box::new(writer))
}