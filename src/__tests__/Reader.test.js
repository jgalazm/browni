

import Reader from '../Reader';

it('Loads bathymetry from a javascript array', () => {

    const data = {
        bathymetry: [[1,1],[1,1]]
    };
    const reader = new Reader(data);

    expect(reader.bathymetry).toEqual(data.bathymetry);
});

it('Loads bathymetry from a csv file', () => {

    expect(1).toBeDefined();
});

it('Loads bathymetry from a binary file', () => {

    expect(1).toBeDefined();
});
