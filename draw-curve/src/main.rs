use std::str::FromStr;

use clap::{Parser, ArgEnum};
use plotters::prelude::*;

use css_color_parser::{Color as CssColor, ColorParseError};

use insides::*;

const MAX_DEPTH: usize = 16;

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ArgEnum)]
enum Curve {
    Morton,
    Hilbert,
}

#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct LocalColor(pub u8, pub u8, pub u8);

impl FromStr for LocalColor {
    type Err = ColorParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let c = s.parse::<CssColor>()?;
        Ok(LocalColor(c.r, c.g, c.b))
    }
}

impl Into<RGBColor> for LocalColor {
    fn into(self) -> RGBColor {
        RGBColor(self.0, self.1, self.2)
    }
}

fn lerp(a: f32, b: f32, t: f32) -> f32 {
    (b - a) * t + a
}

fn lerp_colour(a: LocalColor, b: LocalColor, t: f32) -> LocalColor {
    LocalColor(lerp(a.0 as f32, b.0 as f32, t) as u8, lerp(a.1 as f32, b.1 as f32, t) as u8, lerp(a.2 as f32, b.2 as f32, t) as u8)
}

/// Search for a pattern in a file and display the lines that contain it.
#[derive(Parser)]
#[clap(about = "Produces space filling curve images")]
struct Args {
    #[clap(arg_enum)]
    curve: Curve,

    #[clap(short, long, help = "Curve depth iterations", default_value_t = 1)]
    depth: usize,

    #[clap(short, long, help = "Output image file", default_value = "./curve.png", parse(from_os_str))]
    outfile: std::path::PathBuf,

    #[clap(short = 's', long, help = "Output image edge size in pixels", default_value = "512")]
    outsize: usize,

    #[clap(short = 'r', long, help = "Gif framerate (when outfile extension is .gif)", default_value = "15")]
    framerate: usize,

    #[clap(short = 'f', long, help = "Fade line colour towards end of curve")]
    fade: bool,

    #[clap(short = 'p', long, help = "Opens the image in your default viewer")]
    present: bool,

    #[clap(short = 'a', long, help = "CSS compatible colour A (line colour)", default_value = "rgb(155, 180, 190)")]
    cola: LocalColor,

    #[clap(short = 'b', long, help = "CSS compatible colour B (background colour)", default_value = "rgb(45, 65, 75)")]
    colb: LocalColor,
}

fn draw_curve<T, B>(backend: &mut B, args: &Args) where T: SpaceFillingCurve<2, Coord = u8, Index = u16>, B: DrawingBackend {
    let curve_edge_width = 1 << args.depth;
    let curve_len = curve_edge_width * curve_edge_width;

    let cell_size = args.outsize as f32 / curve_edge_width as f32;
    let half_cell_size = cell_size * 0.5;

    let cell_center = |c: [T::Coord; 2]| ((c[0] as f32 * cell_size + half_cell_size) as i32, args.outsize as i32 - (c[1] as f32 * cell_size + half_cell_size) as i32);

    let blob_radius = half_cell_size * 0.2;
    let line_width = half_cell_size * 0.15;
    let mut line_col: RGBColor = args.cola.into();

    for i in 0..curve_len {
        if args.fade {
            let t = i as f32 / curve_len as f32;
            line_col = lerp_colour(args.cola, args.colb, t).into();
        }
        let sfc = T::from_index(i as T::Index);
        let coords = sfc.coords();
        let b = cell_center(coords);
        let style = ShapeStyle::from(line_col).stroke_width(line_width as u32);
        backend.draw_circle(b, blob_radius as u32, &style, true).expect("Failed to draw circle");
        if i > 0 {
            let a = cell_center(T::from_index(i as T::Index - 1).coords());
            backend.draw_line(a, b, &style).expect("Failed to draw line");
        }
        backend.present().expect("Failed to present frame");
    }
}

fn draw<B>(mut backend: B, args: &Args) where B: DrawingBackend {
    let bg_col: RGBColor = args.colb.into();
    backend.draw_rect(
        (0, 0),
        (args.outsize as i32 - 1, args.outsize as i32 - 1),
        &ShapeStyle::from(bg_col),
        true,
    ).expect("Failed to draw rect");
    backend.present().expect("Failed to present frame");

    match args.curve {
        Curve::Morton => draw_curve::<Morton::<Expand<u8, 2>, 2>, _>(&mut backend, args),
        Curve::Hilbert => draw_curve::<Hilbert::<Expand<u8, 2>, 2>, _>(&mut backend, args),
    };
}

fn main() {
    let args = Args::parse();
    assert!(args.depth >= 1 && args.depth <= MAX_DEPTH, "Parameter 'depth' must have a value between 1 and {MAX_DEPTH} (inclusive)");
    assert!(args.outsize >= 16 && args.outsize <= 4096, "Parameter 'outsize' must have a value between 16 and 4096 pixels (inclusive)");

    if args.outfile.extension().filter(|&e| e == "gif").is_some() {
        let frame_delay = (1000.0 / args.framerate as f32) as u32;
        draw(BitMapBackend::gif(&args.outfile, (args.outsize as u32, args.outsize as u32), frame_delay).expect("Failed to create gif backend"), &args);
    } else {
        draw(BitMapBackend::new(&args.outfile, (args.outsize as u32, args.outsize as u32)), &args);
    }

    if args.present {
        let outfile = args.outfile.canonicalize().expect("Failed to canonicalise outfile");
        opener::open(outfile.clone()).expect(format!("Failed to open file: {}", outfile.to_str().unwrap()).as_str());
    }
}
