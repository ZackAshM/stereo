import stereo.flightpath as flightpath


def try_flight(version: int) -> None:
    """
    Test that we can load ANITA flight data and access all variables.
    """

    # load the flight
    flight = flightpath.load_flight(version)

    # get the number of elements in time
    N = flight.realTime.size

    # check that all the required variables are there
    assert flight.realTime.size == N
    assert flight.altitude.size == N
    assert flight.latitude.size == N
    assert flight.longitude.size == N
    assert flight.heading.size == N

    # and if we are A3/A4, try for the others
    if version in [3, 4]:
        assert flight.pitch.size == N
        assert flight.roll.size == N

    # and we are done


def test_anita4() -> None:
    """
    Test loading ANITA4 flight data.
    """
    try_flight(4)
