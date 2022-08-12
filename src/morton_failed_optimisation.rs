
        let mut index = Self::Index::zero();

        // Calculate how many bits are available
        let coord_bits = DM::UNDILATED_BITS;
        let total_coord_bits = DM::UNDILATED_BITS * D;
        let page_bits = (core::mem::size_of::<usize>() * 8) / D;

        // Iterate all bits in all coords
        let mut coord_cursor = 0;
        while coord_cursor < total_coord_bits {
            let mut page = 0;
            let mut page_cursor = 0;
            let coord_cursor_reset_point = coord_cursor;

            while coord_cursor < total_coord_bits && page_cursor < page_bits {
                // Work out which coordinate we're in and how far in
                let mapped_coord = coord_cursor / coord_bits;
                let mapped_coord_cursor = coord_cursor % coord_bits;

                // Calculate number of bits to copy
                let coord_bits_remaining = coord_bits - mapped_coord_cursor;
                let page_bits_remaining = page_bits - page_cursor;
                let bits_to_copy = if coord_bits_remaining < page_bits_remaining {
                    coord_bits_remaining
                } else {
                    page_bits_remaining
                };
                let lower_undilated_mask = (1 << bits_to_copy) - 1;
                let lower_dilated_mask = build_dilated_mask(bits_to_copy, D);

                // Map bits to page
                page |= (coords[mapped_coord].shr(mapped_coord_cursor).to_usize()
                    & lower_undilated_mask)
                    << page_cursor;

                // Move to next coord or page
                coord_cursor += bits_to_copy;
                page_cursor += bits_to_copy;
            }
            page = dilate_usize::<D>(page);

            let mut page_cursor = 0;
            coord_cursor = coord_cursor_reset_point;

            while coord_cursor < total_coord_bits && page_cursor < page_bits {
                // Work out which coordinate we're in and how far in
                let mapped_coord = coord_cursor / coord_bits;
                let mapped_coord_cursor = coord_cursor % coord_bits;

                // Calculate number of bits to copy
                let coord_bits_remaining = coord_bits - mapped_coord_cursor;
                let page_bits_remaining = page_bits - page_cursor;
                let bits_to_copy = if coord_bits_remaining < page_bits_remaining {
                    coord_bits_remaining
                } else {
                    page_bits_remaining
                };
                let lower_undilated_mask = (1 << bits_to_copy) - 1;
                let lower_dilated_mask = build_dilated_mask(bits_to_copy, D);

                // Map page to index
                index = index.bit_or(
                    Self::Index::from_usize(
                        ((page >> (page_cursor * D)) & lower_dilated_mask) << mapped_coord,
                    )
                    .shl(mapped_coord_cursor * D),
                );

                // Move to next coord or page
                coord_cursor += bits_to_copy;
                page_cursor += bits_to_copy;
            }
        }
        Self(index)

        //        Self(
        //            coords
        //                .into_iter()
        //                .enumerate()
        //                .fold(Self::Index::zero(), |v, (i, c)| {
        //                    v.bit_or(DM::dilate(c).value().shl(i))
        //                }),
        //        )
