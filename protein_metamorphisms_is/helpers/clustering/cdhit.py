def calculate_cdhit_word_length(identity_threshold, logger):
    """
        Determines the word length for CD-HIT based on the identity threshold.

        Parameters
        ----------
        identity_threshold : float
            The sequence identity threshold (between 0 and 1).

        Returns
        -------
        int
            The corresponding word length for CD-HIT clustering.

        Notes
        -----
        - If `identity_threshold >= 0.7`, returns 5.
        - If `0.6 <= identity_threshold < 0.7`, returns 4.
        - If `0.5 <= identity_threshold < 0.6`, returns 3.
        - For all other cases, returns 5 by default.
        """
    logger.info(f"Calculating CD-HIT word length for {identity_threshold}.")

    if identity_threshold >= 0.7:
        word_length = 5

    elif identity_threshold >= 0.6:
        word_length = 4

    elif identity_threshold >= 0.5:
        word_length = 3
    else:
        word_length = 5  # Default case

    logger.info(f"For a threshold of {identity_threshold}, the calculated word length is: {word_length}")
    return word_length
